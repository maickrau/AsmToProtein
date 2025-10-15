#!/usr/bin/env python

import os
import sys
import pathlib
import subprocess
import sqlite3
import shutil
import datetime
import itertools
import SequenceReader
import Gff3Parser
import HandleAssembly
import Util

def isoform_sort_order(isoform):
	"""
	Orders isoforms so ref is first, rest ordered by A, B, C, ... ZZ, AA, AB, and novel isoforms afterwards ordered by novel_A, novel_B ...

	Args:
		isoform: Isoform name

	Returns:
		Integer describing isoform order
	"""
	if isoform == "ref": return 0
	is_novel = False
	if isoform[0:6] == "novel_":
		is_novel = True
		isoform = isoform[6:]
	result = 1
	for c in isoform[::-1]:
		result *= 26
		result += ord(c) - ord('A') + 1
	if is_novel:
		result *= 10000 # hope there aren't more than 10000 isoforms per transcript lol
	return result

def alleleset_sort_order(alleleset):
	"""
	Orders an allele set so missing is first, followed by elementwise comparisons
	"""
	if len(alleleset) == 0:
		return (0)
	result = []
	for c in alleleset:
		result.append(isoform_sort_order(c) + 1)
	result = tuple(result)
	return result

def get_isoform_name(index):
	result = ""
	while True:
		result += "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[int(index % 26)]
		if index < 26: break
		index = int(index/26)
	return "".join(result[::-1])

def remove_temporary_groups(database_path, groups):
	"""
	Removes temporary groups.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		groups: List of group names
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		cursor.execute(f"DELETE FROM SampleGroup WHERE GroupName IN ({",".join("?" for group in groups)})", groups)

def add_temporary_groups(database_path, groups):
	"""
	Creates temporary groups which contain specified samples

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		groups: Dictionary of sample -> group

	Returns:
		Dictionary of name mappings created_temporary_group_name -> original_group_name
	"""
	name_mapping = {}
	for _, group in groups.items():
		if group not in name_mapping: name_mapping[group] = "tmp_" + Util.make_random_prefix()
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		for original_name in name_mapping:
			# make sure temp names are unique and unused, generate new ones until so
			while True:
				temp_name = name_mapping[original_name]
				found = cursor.execute(f"SELECT COUNT(*) FROM SampleGroup WHERE GroupName=?", (temp_name,)).fetchone()
				if found[0] == 0:
					break
				name_mapping[original_name] = "tmp_" + Util.make_random_prefix()
				continue
		additions = []
		sample_db_id = {}
		for sample, db_id in cursor.execute("SELECT Name, Id FROM Sample").fetchall():
			sample_db_id[sample] = db_id
		for sample, group in groups.items():
			if sample not in sample_db_id: raise RuntimeError(f"Sample {sample} not found. Note that sample names are case sensitive")
			additions.append((sample_db_id[sample], name_mapping[group]))
		cursor.executemany("INSERT INTO SampleGroup (SampleId, GroupName) VALUES (?, ?)", additions)
	result = {}
	for real_name, temporary_name in name_mapping.items():
		result[temporary_name] = real_name
	return result

def check_if_haplotypes_are_fine(database_path):
	"""
	Checks if haplotypes are as expected with two haplotypes per each sample, either 1&2 or mat&pat

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of errors. If no errors, empty list
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		sample_haplotypes = {}
		result = []
		for sample, haplotype in cursor.execute("SELECT Sample.Name, Haplotype.Haplotype FROM Sample INNER JOIN Haplotype ON Haplotype.SampleId=Sample.Id"):
			if sample not in sample_haplotypes: sample_haplotypes[sample] = []
			sample_haplotypes[sample].append(haplotype)
		for sample in sample_haplotypes:
			if sample == "reference":
				if len(sample_haplotypes[sample]) != 1:
					result.append(f"Reference should have exactly one haplotype, now has {len(sample_haplotypes[sample])}")
				elif sample_haplotypes[sample][0] != "reference":
					result.append(f"Reference should have haplotype with name \"reference\", now has \"{sample_haplotypes[sample][0]}\"")
			else:
				sample_haplotypes[sample].sort()
				if len(sample_haplotypes[sample]) != 2:
					result.append(f"Sample \"{sample}\" has an unexpected number of haplotypes: {len(sample_haplotypes[sample])} (haplotypes are {", ".join("\"" + hap + "\"" for hap in sample_haplotypes[sample])}), expected 2 haplotypes")
				else:
					valid = False
					if sample_haplotypes[sample][0] == "1" and sample_haplotypes[sample][1] == "2":
						valid = True
					if sample_haplotypes[sample][0] == "mat" and sample_haplotypes[sample][1] == "pat":
						valid = True
					if not valid:
						result.append(f"Sample \"{sample}\" has unexpected haplotypes, expected either \"1 & 2\" or \"mat & pat\", sample has \"{" & ".join(sample_haplotypes[sample])}\"")
		return result

def basic_stats(database_path):
	"""
	Print a bunch of basic stats about the database

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		sample_count = cursor.execute("SELECT COUNT(*) FROM Sample").fetchall()[0][0]
		print(f"Number of samples including reference: {sample_count}")
		print(f"Number of samples excluding reference: {sample_count-1}")
		haplotype_count = cursor.execute("SELECT COUNT(*) FROM Haplotype").fetchall()[0][0]
		print(f"Number of haplotypes including reference: {haplotype_count}")
		print(f"Number of haplotypes excluding reference: {haplotype_count-1}")
		gene_count = cursor.execute("SELECT COUNT(*) FROM Gene").fetchall()[0][0]
		print(f"Number of genes: {gene_count}")
		transcript_count = cursor.execute("SELECT COUNT(*) FROM Transcript").fetchall()[0][0]
		print(f"Number of transcripts: {transcript_count}")
		isoform_count = cursor.execute("SELECT COUNT(*) FROM Isoform").fetchall()[0][0]
		print(f"Number of distinct isoforms: {isoform_count}")
		sample_allele_count = cursor.execute("SELECT COUNT(*) FROM SampleAllele").fetchall()[0][0]
		print(f"Number of alleles in all samples: {sample_allele_count}")

def get_all_transcripts(database_path):
	"""
	Gets all transcript names

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of transcript names
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		transcripts = list(x[0] for x in cursor.execute("SELECT Name FROM Transcript ORDER BY Name", ()).fetchall())
		return transcripts

def add_isoforms_to_fasta(file_path, novel_isoforms):
	"""
	Adds isoform sequences to isoform fasta

	Args:
		file_path: Path to isoform fasta
		novel_isoforms: Set of novel isoforms in format (transcript_id, isoform_sequence)

	Returns:
		List of inserted isoforms in format [(transcript_id, sequence_id)]
	"""
	next_free_sequence_id = 0
	for name, _ in SequenceReader.stream_sequences(file_path):
		isoform_id = int(name)
		next_free_sequence_id = max(next_free_sequence_id, isoform_id)
	next_free_sequence_id += 1
	result = []
	with open(file_path, "a") as f:
		for transcript_id, sequence in novel_isoforms:
			f.write(">" + str(next_free_sequence_id) + "\n")
			f.write(sequence + "\n")
			result.append((transcript_id, next_free_sequence_id))
			next_free_sequence_id += 1
	return result

def get_isoform_sequences_from_file(file_path):
	"""
	Gets all isoform sequences

	Args:
		file_path: Path to isoform fasta

	Returns:
		Dictionary of id -> sequence
	"""
	result = {}
	for name, sequence in SequenceReader.stream_sequences(file_path):
		result[int(name)] = sequence
	return result

def get_isoform_sequences_per_transcript(base_path):
	"""
	Gets all isoform sequences stratified by transcript

	Args:
		base_path: Path to database folder. Type should be pathlib.Path

	Returns:
		Dictionary of transcript -> [list of (isoform_sequence, frequency)]
	"""

	isoform_sequences = get_isoform_sequences_from_file(str(base_path / "isoforms.fa"))

	with sqlite3.connect(str(base_path / "sample_info.db")) as connection:
		cursor = connection.cursor()
		result = {}
		for transcript, sequenceID, count in cursor.execute("SELECT Transcript.Name, Isoform.SequenceId, COUNT(*) FROM SampleAllele INNER JOIN Isoform ON SampleAllele.IsoformId=Isoform.Id INNER JOIN Transcript ON Isoform.TranscriptId=Transcript.Id GROUP BY Transcript.Name, Isoform.SequenceID"):
			if transcript not in result: result[transcript] = []
			result[transcript].append((isoform_sequences[sequenceID], count))
		return result

def add_sample_to_group(database_path, sample_name, sample_group):
	"""
	Adds a sample to a group.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		sample_name: Name of sample
		sample_group: Name of group
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		sample_current_groups = set(cursor.execute("SELECT SampleGroup.GroupName FROM SampleGroup INNER JOIN Sample ON SampleGroup.SampleId=Sample.Id WHERE Sample.Name=?", (sample_name,)).fetchall())
		if (sample_group,) in sample_current_groups:
			# already in group, no need to do anything
			return
		sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
		assert sample_id is not None
		sample_id = sample_id[0]
		cursor.execute("INSERT INTO SampleGroup (SampleId, GroupName) VALUES (?, ?)", (sample_id, sample_group))
		connection.commit()

def remove_sample_from_group(database_path, sample_name, sample_group):
	"""
	Removes a sample from a group.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		sample_name: Name of sample
		sample_group: Name of group
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
		assert sample_id is not None
		sample_id = sample_id[0]
		cursor.execute("DELETE FROM SampleGroup WHERE SampleId=? AND GroupName=?", (sample_id, sample_group))
		connection.commit()

def get_transcript_gene_chromosome_info(database_path):
	"""
	Gets transcript's genes and reference chromosomes for all transcripts

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		Dict of transcript_id -> (gene_common_name, reference_chromosome)
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		result = {}
		for transcript_name, gene_common_name, ref_chrom in cursor.execute("SELECT Transcript.Name, Gene.CommonName, Gene.ReferenceChromosome FROM Transcript INNER JOIN Gene ON Transcript.GeneId=Gene.Id", ()):
			result[transcript_name] = (gene_common_name, ref_chrom)
		return result

def get_isoforms_of_transcript(base_path, transcript):
	"""
	Gets all isoforms of a single transcript

	Args:
		base_path: Path to database folder. Type should be pathlib.Path
		transcript: Name of transcript

	Returns:
		List of [(isoform_name, isoform_sequence, isoform_copy_count)]
	"""
	
	isoform_sequences = get_isoform_sequences_from_file(str(base_path / "isoforms.fa"))

	with sqlite3.connect(str(base_path / "sample_info.db")) as connection:
		cursor = connection.cursor()
		result = []
		for name, sequenceID, copycount in cursor.execute("SELECT Isoform.Name, Isoform.SequenceID, COUNT(*) FROM Isoform INNER JOIN Transcript ON Isoform.TranscriptId=Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId=Isoform.Id WHERE Transcript.Name=? GROUP BY Isoform.Name, Isoform.SequenceID", (transcript,)):
			result.append((name, isoform_sequences[sequenceID], copycount))
		result.sort()
		return result

def get_isoforms_of_all_transcripts(base_path):
	"""
	Gets all isoforms of all transcripts

	Args:
		base_path: Path to database folder. Type should be pathlib.Path

	Returns:
		List of [(transcript, list of [(isoform_name, isoform_sequence, isoform_copy_count)])]
	"""

	isoform_sequences = get_isoform_sequences_from_file(str(base_path / "isoforms.fa"))

	with sqlite3.connect(str(base_path / "sample_info.db")) as connection:
		cursor = connection.cursor()
		result_per_transcript = {}
		isoform_coverage = {}
		for isoformid, coverage in cursor.execute("SELECT IsoformId, COUNT(*) FROM SampleAllele GROUP BY IsoformId", ()):
			isoform_coverage[isoformid] = coverage
		for transcript, name, sequenceID, isoformid in cursor.execute("SELECT Transcript.Name, Isoform.Name, Isoform.SequenceID, Isoform.Id FROM Isoform INNER JOIN Transcript ON Isoform.TranscriptId=Transcript.Id", ()):
			if transcript not in result_per_transcript: result_per_transcript[transcript] = []
			result_per_transcript[transcript].append((name, isoform_sequences[sequenceID], isoform_coverage[isoformid]))
		result = []
		for transcript, isoforms in result_per_transcript.items():
			isoforms.sort()
			result.append((transcript, isoforms))
		result.sort()
		return result

def get_group_sample_counts(database_path, groups):
	"""
	Gets counts of samples in groups.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		groups: List of group names

	Returns:
		Dictionary of (group -> count)
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		result = {}
		for group, count in cursor.execute(f"SELECT GroupName, COUNT(*) FROM SampleGroup WHERE GroupName IN ({",".join("?" for group in groups)}) GROUP BY GroupName", groups):
			result[group] = count
		return result

def get_all_allelesets_per_haplotype(database_path):
	"""
	Gets allele sets of all transcripts of all samples per haplotype.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of (transcript, list of (sample_name, sample_haplotype, alleleset))
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		allelesets_per_sample = {}
		all_samples = set()
		for sample, haplotype in cursor.execute("SELECT Sample.Name, Haplotype.Haplotype FROM Sample INNER JOIN Haplotype ON Sample.Id=Haplotype.SampleId"):
			all_samples.add((sample, haplotype))
		all_samples = list(all_samples)
		all_samples.sort()
		for transcript, sample, haplotype, isoform in cursor.execute(f"SELECT Transcript.Name, Sample.Name, Haplotype.Haplotype, Isoform.Name FROM Transcript INNER JOIN Isoform ON Isoform.TranscriptId = Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN Sample ON Haplotype.SampleId=Sample.Id", ()).fetchall():
			if transcript not in allelesets_per_sample: allelesets_per_sample[transcript] = {}
			if (sample, haplotype) not in allelesets_per_sample[transcript]: allelesets_per_sample[transcript][(sample, haplotype)] = []
			allelesets_per_sample[transcript][(sample, haplotype)].append(isoform)
		for transcript, samples in allelesets_per_sample.items():
			for (sample, haplotype), alleles in samples.items():
				alleleset = alleles
				alleleset.sort(key=lambda x: isoform_sort_order(x))
				alleleset = tuple(alleleset)
				allelesets_per_sample[transcript][(sample, haplotype)] = alleleset
		result = []
		for transcript in allelesets_per_sample:
			sample_lines = []
			for (sample, haplotype) in all_samples:
				alleleset = allelesets_per_sample[transcript].get((sample, haplotype), ())
				sample_lines.append((sample, haplotype, Util.get_alleleset_name(alleleset)))
			result.append((transcript, sample_lines))
		result.sort()
		return result

def get_allelesets_of_one_transcript(database_path, transcript):
	"""
	Gets allele sets of one transcripts of all samples.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of (transcript, list of (sample_name, alleleset))
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		allelesets_per_sample = {}
		all_samples = set()
		for sample, in cursor.execute("SELECT Sample.Name FROM Sample"):
			all_samples.add(sample)
		all_samples = list(all_samples)
		all_samples.sort()
		for transcript, sample, isoform in cursor.execute(f"SELECT Transcript.Name, Sample.Name, Isoform.Name FROM Transcript INNER JOIN Isoform ON Isoform.TranscriptId = Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN Sample ON Haplotype.SampleId=Sample.Id WHERE Transcript.Name=?", (transcript,)).fetchall():
			if transcript not in allelesets_per_sample: allelesets_per_sample[transcript] = {}
			if sample not in allelesets_per_sample[transcript]: allelesets_per_sample[transcript][sample] = []
			allelesets_per_sample[transcript][sample].append(isoform)
		for transcript, samples in allelesets_per_sample.items():
			for sample, alleles in samples.items():
				alleleset = alleles
				alleleset.sort(key=lambda x: isoform_sort_order(x))
				alleleset = tuple(alleleset)
				allelesets_per_sample[transcript][sample] = alleleset
		result = []
		for transcript in allelesets_per_sample:
			sample_lines = []
			for sample in all_samples:
				alleleset = allelesets_per_sample[transcript].get(sample, ())
				sample_lines.append((sample, Util.get_alleleset_name(alleleset)))
			result.append((transcript, sample_lines))
		result.sort()
		return result

def get_all_allelesets(database_path):
	"""
	Gets allele sets of all transcripts of all samples.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of (transcript, list of (sample_name, alleleset))
	"""
	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		allelesets_per_sample = {}
		all_samples = set()
		for sample, in cursor.execute("SELECT Sample.Name FROM Sample"):
			all_samples.add(sample)
		all_samples = list(all_samples)
		all_samples.sort()
		for transcript, sample, isoform in cursor.execute(f"SELECT Transcript.Name, Sample.Name, Isoform.Name FROM Transcript INNER JOIN Isoform ON Isoform.TranscriptId = Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN Sample ON Haplotype.SampleId=Sample.Id", ()).fetchall():
			if transcript not in allelesets_per_sample: allelesets_per_sample[transcript] = {}
			if sample not in allelesets_per_sample[transcript]: allelesets_per_sample[transcript][sample] = []
			allelesets_per_sample[transcript][sample].append(isoform)
		for transcript, samples in allelesets_per_sample.items():
			for sample, alleles in samples.items():
				alleleset = alleles
				alleleset.sort(key=lambda x: isoform_sort_order(x))
				alleleset = tuple(alleleset)
				allelesets_per_sample[transcript][sample] = alleleset
		result = []
		for transcript in allelesets_per_sample:
			sample_lines = []
			for sample in all_samples:
				alleleset = allelesets_per_sample[transcript].get(sample, ())
				sample_lines.append((sample, Util.get_alleleset_name(alleleset)))
			result.append((transcript, sample_lines))
		result.sort()
		return result

def get_samples(database_path):
	"""
	Lists all samples and their haplotypes.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of (sample_name, list of (haplotype_name))
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		result_per_sample = {}
		for sample, haplotype in cursor.execute("SELECT Sample.Name, Haplotype.Haplotype FROM Haplotype INNER JOIN Sample ON Haplotype.SampleId=Sample.Id"):
			if sample not in result_per_sample: result_per_sample[sample] = set()
			result_per_sample[sample].add(haplotype)
		result = []
		for sample, haplotypes in result_per_sample.items():
			sorted_haps = list(haplotypes)
			sorted_haps.sort()
			result.append((sample, sorted_haps))
		result.sort()
		return result

def get_sample_groups(database_path):
	"""
	Lists all samples and their (possibly zero) groups.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path

	Returns:
		List of (sample_name, list of (group_name))
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		result_per_sample = {}
		for sample, group in cursor.execute("SELECT Sample.Name, SampleGroup.GroupName FROM Sample LEFT JOIN SampleGroup ON SampleGroup.SampleId=Sample.Id"):
			if sample not in result_per_sample: result_per_sample[sample] = set()
			if group:
				result_per_sample[sample].add(group)
		result = []
		for sample, groups in result_per_sample.items():
			sorted_groups = list(groups)
			sorted_groups.sort()
			result.append((sample, sorted_groups))
		result.sort()
		return result

def get_all_transcripts_alleleset_contingency_table_by_group(database_path, groups):
	"""
	Gets contingency table of allele sets with samples divided into first_group and second_group.
	If samples are shared between groups, throws error.
	All transcripts included.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		groups: List of group names

	Returns:
		List of (transcript, list of (alleleset, [count_in_group_n]))
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		samples_per_group = {}
		groups_per_sample = {}
		# check for any overlapping samples
		for sample, group in cursor.execute(f"SELECT SampleId, GroupName FROM SampleGroup WHERE GroupName IN ({",".join("?" for group in groups)})", tuple(groups)).fetchall():
			if group not in samples_per_group: samples_per_group[group] = set()
			samples_per_group[group].add(sample)
			if sample not in groups_per_sample: groups_per_sample[sample] = set()
			groups_per_sample[sample].add(group)
		for sample in groups_per_sample:
			if len(groups_per_sample[sample]) >= 2:
				raise RuntimeError(f"Samples are shared between groups {", ".join(groups_per_sample[sample])}")
		allelesets_per_sample = {}
		for transcript, in cursor.execute(f"SELECT Transcript.Name FROM Transcript", ()).fetchall():
			if transcript not in allelesets_per_sample: allelesets_per_sample[transcript] = {}
		for transcript, sample, isoform in cursor.execute(f"SELECT Transcript.Name, Haplotype.SampleId, Isoform.Name FROM Transcript INNER JOIN Isoform ON Isoform.TranscriptId = Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN SampleGroup ON Haplotype.SampleId=SampleGroup.SampleId WHERE SampleGroup.GroupName IN ({",".join("?" for group in groups)})", tuple(groups)).fetchall():
			assert transcript in allelesets_per_sample
			if sample not in allelesets_per_sample[transcript]: allelesets_per_sample[transcript][sample] = []
			allelesets_per_sample[transcript][sample].append(isoform)
		for transcript, samples in allelesets_per_sample.items():
			for sample, alleles in samples.items():
				alleleset = alleles
				alleleset.sort(key=lambda x: isoform_sort_order(x))
				alleleset = tuple(alleleset)
				allelesets_per_sample[transcript][sample] = alleleset
		counts = {}
		for transcript in allelesets_per_sample:
			counts[transcript] = {}
			for group, samples in samples_per_group.items():
				for sample in samples:
					if sample not in allelesets_per_sample[transcript]:
						allelesets_per_sample[transcript][sample] = () # empty tuple
					if allelesets_per_sample[transcript][sample] not in counts[transcript]: counts[transcript][allelesets_per_sample[transcript][sample]] = [0] * len(groups)
		for transcript in allelesets_per_sample:
			for i in range(0, len(groups)):
				group = groups[i]
				for sample in samples_per_group[group]:
					counts[transcript][allelesets_per_sample[transcript][sample]][i] += 1
		result = []
		for transcript, countinfo in counts.items():
			lines_here = []
			for alleleset, groupcounts in countinfo.items():
				lines_here.append((Util.get_alleleset_name(alleleset), groupcounts))
			lines_here.sort(key=lambda x: alleleset_sort_order(x[0]))
			result.append((transcript, lines_here))
		result.sort()
		return result

def get_alleleset_contingency_table_by_group(database_path, transcript_name, groups):
	"""
	Gets contingency table of allele sets with samples divided into first_group and second_group.
	If samples are shared between groups, throws error

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		transcript_name: Name of transcript
		groups: List of group names

	Returns:
		List of (alleleset, [count_in_group_n])
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		samples_per_group = {}
		groups_per_sample = {}
		# check for any overlapping samples
		for sample, group in cursor.execute(f"SELECT SampleId, GroupName FROM SampleGroup WHERE GroupName IN ({",".join("?" for group in groups)})", tuple(groups)).fetchall():
			if group not in samples_per_group: samples_per_group[group] = set()
			samples_per_group[group].add(sample)
			if sample not in groups_per_sample: groups_per_sample[sample] = set()
			groups_per_sample[sample].add(group)
		for sample in groups_per_sample:
			if len(groups_per_sample[sample]) >= 2:
				raise RuntimeError(f"Samples are shared between groups {", ".join(groups_per_sample[sample])}")
		allelesets_per_sample = {}
		for sample, isoform in cursor.execute(f"SELECT Haplotype.SampleId, Isoform.Name FROM Transcript INNER JOIN Isoform ON Isoform.TranscriptId = Transcript.Id INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN SampleGroup ON Haplotype.SampleId=SampleGroup.SampleId WHERE Transcript.Name = ? AND SampleGroup.GroupName IN ({",".join("?" for group in groups)})", tuple([transcript_name] + groups)).fetchall():
			if sample not in allelesets_per_sample: allelesets_per_sample[sample] = []
			allelesets_per_sample[sample].append(isoform)
		for sample in allelesets_per_sample:
			alleleset = allelesets_per_sample[sample]
			alleleset.sort(key=lambda x: isoform_sort_order(x))
			alleleset = tuple(alleleset)
			allelesets_per_sample[sample] = alleleset
		counts = {}
		for group, samples in samples_per_group.items():
			for sample in samples:
				if sample not in allelesets_per_sample:
					allelesets_per_sample[sample] = () # empty tuple
				if allelesets_per_sample[sample] not in counts: counts[allelesets_per_sample[sample]] = [0] * len(groups)
		for i in range(0, len(groups)):
			group = groups[i]
			for sample in samples_per_group[group]:
				counts[allelesets_per_sample[sample]][i] += 1
		result = []
		for alleleset, groupcounts in counts.items():
			result.append((Util.get_alleleset_name(alleleset), groupcounts))
		result.sort(key=lambda x: alleleset_sort_order(x[0]))
		return result

def rename_isoforms_by_coverage(database_path):
	"""
	Renames the isoforms in the database according to their coverage. Each transcript's isoforms have their own namespaces, names are ordered A, B, C, ... Z, AA, AB, ... ZZ, AAA, AAB, ...
	Reference isoform gets special name "ref", unless reference has multiple isoforms in which case they have no special name

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		isoform_coverages_per_transcript = {}
		reference_isoform = {}
		isoform_belongs_to_transcript = {}
		print(f"{datetime.datetime.now().astimezone()}: Get isoform to transcript mapping", file=sys.stderr)
		for transcript, isoform in cursor.execute("SELECT TranscriptId, Id FROM Isoform"):
			isoform_belongs_to_transcript[isoform] = transcript
		print(f"{datetime.datetime.now().astimezone()}: Get isoform coverages", file=sys.stderr)
		for isoform, coverage in cursor.execute("SELECT IsoformId, COUNT(*) FROM SampleAllele GROUP BY IsoformId").fetchall():
			transcript = isoform_belongs_to_transcript[isoform]
			if transcript not in isoform_coverages_per_transcript: isoform_coverages_per_transcript[transcript] = []
			isoform_coverages_per_transcript[transcript].append((coverage, isoform))
		print(f"{datetime.datetime.now().astimezone()}: Get reference isoforms", file=sys.stderr)
		for transcript, isoform in cursor.execute("SELECT DISTINCT TranscriptId, Isoform.Id FROM Isoform INNER JOIN SampleAllele ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN Sample ON Haplotype.SampleId = Sample.Id WHERE Sample.Name='reference'").fetchall():
			if transcript not in reference_isoform:
				reference_isoform[transcript] = isoform
			else:
				reference_isoform[transcript] = None
		print(f"{datetime.datetime.now().astimezone()}: Update reference isoform names", file=sys.stderr)
		for transcript, isoform in reference_isoform.items():
			if isoform is None: continue
			cursor.execute("UPDATE Isoform SET Name=? WHERE Id=?", ("ref", isoform))
		print(f"{datetime.datetime.now().astimezone()}: Update alt isoform names", file=sys.stderr)
		for transcript, coveragelist in isoform_coverages_per_transcript.items():
			coveragelist.sort(key=lambda x: (-x[0], x[1]))
			index = 0
			for coverage, isoform in coveragelist:
				if isoform == reference_isoform[transcript]: continue
				name = get_isoform_name(index)
				cursor.execute("UPDATE Isoform SET Name=? WHERE Id=?", (name, isoform))
				index += 1
		print(f"{datetime.datetime.now().astimezone()}: Commit", file=sys.stderr)
		connection.commit()
	print(f"{datetime.datetime.now().astimezone()}: Done", file=sys.stderr)

def initialize_sqlite_db_schema(database_path):
	"""
	Initializes the sample annotation sql database schema

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
	"""

	make_tables_command = \
	"""
		CREATE TABLE Sample (Id INTEGER PRIMARY KEY NOT NULL, Name TEXT NOT NULL);
		CREATE TABLE SampleGroup (Id INTEGER PRIMARY KEY NOT NULL, SampleId INTEGER NOT NULL, GroupName TEXT NOT NULL, FOREIGN KEY (SampleId) REFERENCES Sample(Id));
		CREATE TABLE Haplotype (Id INTEGER PRIMARY KEY NOT NULL, SampleId INTEGER NOT NULL, Haplotype TEXT NOT NULL, FOREIGN KEY (SampleId) REFERENCES Sample(Id));
		CREATE TABLE SampleContig (Id INTEGER PRIMARY KEY NOT NULL, HaplotypeId INTEGER NOT NULL, Length INTEGER NOT NULL, Name TEXT NOT NULL, FOREIGN KEY (HaplotypeId) REFERENCES Haplotype(Id));
		CREATE TABLE Gene (Id INTEGER PRIMARY KEY NOT NULL, Name TEXT NOT NULL, CommonName TEXT NOT NULL, ReferenceChromosome TEXT NOT NULL, ReferenceLocationStrand TEXT NOT NULL, ReferenceLocationStart INTEGER NOT NULL, ReferenceLocationEnd INTEGER NOT NULL);
		CREATE TABLE Transcript (Id INTEGER PRIMARY KEY NOT NULL, Name TEXT NOT NULL, GeneId INTEGER NOT NULL, ReferenceChromosome TEXT NOT NULL, ReferenceLocationStrand TEXT NOT NULL, ReferenceLocationStart INTEGER NOT NULL, ReferenceLocationEnd INTEGER NOT NULL, FOREIGN KEY (GeneId) REFERENCES Gene(Id));
		CREATE TABLE Isoform (Id INTEGER PRIMARY KEY NOT NULL, TranscriptId INTEGER NOT NULL, SequenceID INTEGER NOT NULL, Name TEXT, FOREIGN KEY (TranscriptId) REFERENCES Transcript(Id));
		CREATE TABLE SampleAllele (Id INTEGER PRIMARY KEY NOT NULL, HaplotypeId INTEGER NOT NULL, IsoformId INTEGER NOT NULL, SampleContigId INTEGER NOT NULL, SampleLocationStrand TEXT NOT NULL, SampleLocationStart INTEGER NOT NULL, SampleLocationEnd INTEGER NOT NULL, FOREIGN KEY (HaplotypeId) REFERENCES Haplotype(Id), FOREIGN KEY (IsoformId) REFERENCES Isoform(Id), FOREIGN KEY (SampleContigId) REFERENCES SampleContig(Id));
		CREATE TABLE SampleAnnotation (Id INTEGER PRIMARY KEY NOT NULL, SampleId INTEGER NOT NULL, Key TEXT NOT NULL, Value TEXT, FOREIGN KEY (SampleId) REFERENCES Sample(Id));
		CREATE TABLE GeneAnnotation (Id INTEGER PRIMARY KEY NOT NULL, GeneId INTEGER NOT NULL, Key TEXT NOT NULL, Value TEXT, FOREIGN KEY (GeneId) REFERENCES Gene(Id));
		CREATE TABLE TranscriptAnnotation (Id INTEGER PRIMARY KEY NOT NULL, TranscriptId INTEGER NOT NULL, Key TEXT NOT NULL, Value TEXT, FOREIGN KEY (TranscriptId) REFERENCES Transcript(Id));
		CREATE TABLE IsoformAnnotation (Id INTEGER PRIMARY KEY NOT NULL, IsoformId INTEGER NOT NULL, Key TEXT NOT NULL, Value TEXT, FOREIGN KEY (IsoformId) REFERENCES Isoform(Id));
		CREATE INDEX SampleNameIndex ON Sample(Name);
		CREATE INDEX GeneNameIndex ON Gene(Name);
		CREATE INDEX IsoformIndex ON Isoform(TranscriptId, SequenceID);
		CREATE INDEX SampleAlleleIsoformId ON SampleAllele(IsoformId);
		CREATE INDEX SampleAlleleHaplotypeId ON SampleAllele(HaplotypeId);
		CREATE INDEX HaplotypeSampleIdIndex ON Haplotype(SampleId);
		CREATE INDEX IsoformTranscriptId ON Isoform(TranscriptId);
		CREATE INDEX TranscriptNameIndex ON Transcript(Name);
		CREATE INDEX SampleGroupIndex ON SampleGroup(SampleId, GroupName);
		CREATE INDEX SampleGroupIdIndex ON SampleGroup(SampleId);
		CREATE INDEX SampleGroupNameIndex ON SampleGroup(GroupName);
	"""

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		cursor.executescript(make_tables_command)
		connection.commit()

def check_isoforms_have_names(database_folder):
	"""
	Checks that all isoforms have names

	Args:
		database_folder: Path to database folder. Type should be pathlib.Path

	Returns:
		True if all isoforms have names, False otherwise
	"""
	with sqlite3.connect(str(database_folder + "/sample_info.db")) as connection:
		cursor = connection.cursor()
		any_unnamed = cursor.execute("SELECT Id FROM Isoform WHERE Name IS NULL LIMIT 1").fetchall()
		if len(any_unnamed) > 0: return False
	return True

def insert_reference_annotations_to_db(database_path, reference_gff3):
	"""
	Fills the database tables Gene, Transcript with the reference annotation.

	Args:
		database_path: Path to sql database file. Type should be pathlib.Path
		database_path: Path to reference gff3 annotation. Type should be pathlib.Path
	"""

	transcript_genes = Gff3Parser.parse_gff3_gene_names_of_transcripts(reference_gff3)
	(gene_locations, transcript_locations) = Gff3Parser.parse_gene_transcript_locations(reference_gff3)
	
	all_gene_names = set()
	for transcript, (gene_name, gene_id) in transcript_genes.items():
		all_gene_names.add((gene_name, gene_id))

	with sqlite3.connect(str(database_path)) as connection:
		cursor = connection.cursor()
		for (gene_name, gene_id) in all_gene_names:
			cursor.execute("INSERT INTO Gene (Name, CommonName, ReferenceChromosome, ReferenceLocationStrand, ReferenceLocationStart, ReferenceLocationEnd) VALUES (?, ?, ?, ?, ?, ?)", (gene_id, gene_name, gene_locations[gene_id][0], gene_locations[gene_id][1], gene_locations[gene_id][2], gene_locations[gene_id][3]))
		gene_ID_names = {}
		for gene_id, name in cursor.execute("SELECT Id, Name FROM Gene").fetchall():
			assert name not in gene_ID_names
			gene_ID_names[name] = gene_id
		for transcript, (_, gene_id) in transcript_genes.items():
			cursor.execute("INSERT INTO Transcript (Name, GeneId, ReferenceChromosome, ReferenceLocationStrand, ReferenceLocationStart, ReferenceLocationEnd) VALUES (?, ?, ?, ?, ?, ?)", (transcript, gene_ID_names[gene_id], transcript_locations[transcript][0], transcript_locations[transcript][1], transcript_locations[transcript][2], transcript_locations[transcript][3]))
		connection.commit()

def initialize_database(folder, reference_fasta, reference_annotation, liftoff_path, agc_path, num_threads):
	"""
	Creates a database at the specified folder. Runs liftoff on the reference, creates agc read storage, initializes sqlite db.

	Args:
		folder: Target folder for database. Should not exist previously. Type should be pathlib.Path.
		reference_fasta: Fasta file of reference sequence
		reference_annotation: gff3 file of reference annotation
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
		num_threads: Number of threads for liftoff
	"""

	os.makedirs(folder, exist_ok=False)
	refseq_md5sum_check = subprocess.run(["md5sum", reference_fasta], capture_output=True, text=True)
	if refseq_md5sum_check.returncode != 0:
		raise RuntimeError("Could not check reference file md5sum.")
	refannotation_md5sum_check = subprocess.run(["md5sum", reference_annotation], capture_output=True, text=True)
	if refannotation_md5sum_check.returncode != 0:
		raise RuntimeError("Could not check annotation file md5sum.")
	with open(folder / "info.txt", "w") as f:
		print(f"Created on {datetime.datetime.now().astimezone()}", file=f)
		print(f"Reference fasta: {reference_fasta}", file=f)
		print(f"Reference annotation: {reference_annotation}", file=f)
		print(f"Reference fasta md5 checksum: {str(refseq_md5sum_check.stdout)[0:32]}", file=f)
		print(f"Reference annotation md5 checksum: {str(refannotation_md5sum_check.stdout)[0:32]}", file=f)
		liftoff_version = Util.get_liftoff_version(liftoff_path)
		print(f"IsoformCheck version: {Util.Version}", file=f)
		print(f"IsoformCheck db version: {Util.DBVersion}", file=f)
		print(f"Liftoff version: {liftoff_version}", file=f)

	with open(folder / "isoforms.fa", "w") as f:
		# just touch the isoform file to create it
		pass

	tmp_folder = folder / "tmp"
	sample_annotation_folder = folder / "sample_annotations"
	os.makedirs(tmp_folder)
	os.makedirs(sample_annotation_folder)

	ref_fasta = folder / "reference.fa"
	HandleAssembly.prepare_fasta(reference_fasta, ref_fasta)
	print(f"Filtering reference annotation", file=sys.stderr)
	Gff3Parser.filter_gff3_to_things_with_CDS(reference_annotation, str(folder / "reference.gff3"))

	parsed_ref_md5sum_check = subprocess.run(["md5sum", str(folder / "reference.fa")], capture_output=True, text=True)
	if parsed_ref_md5sum_check.returncode != 0:
		raise RuntimeError("Could not check parsed reference fasta md5sum.")
	parsed_ref_hash = str(parsed_ref_md5sum_check.stdout)[0:32]
	filtered_refannotation_md5sum_check = subprocess.run(["md5sum", str(folder / "reference.gff3")], capture_output=True, text=True)
	if filtered_refannotation_md5sum_check.returncode != 0:
		raise RuntimeError("Could not check filtered annotation file md5sum.")
	filtered_refannotation_hash = str(filtered_refannotation_md5sum_check.stdout)[0:32]
	with open(folder / "info.txt", "a") as f:
		print(f"Parsed reference fasta md5 checksum: {parsed_ref_hash}", file=f)
		print(f"Filtered reference annotation md5 checksum: {filtered_refannotation_hash}", file=f)

	liftoff_command = [liftoff_path, "-g", str(folder / "reference.gff3"), "-o", str(tmp_folder / "reference_annotation.gff3"), "-p", str(num_threads), "-sc", "0.95", "-copies", "-dir", str(tmp_folder / "intermediate_files"), "-u", str(tmp_folder / "unmapped_features.txt"), str(ref_fasta), str(ref_fasta)]
	print(f"Running liftoff with command:", file=sys.stderr)
	print(f"{' '.join(liftoff_command)}", file=sys.stderr)
	liftoff_result = subprocess.run(liftoff_command)
	if liftoff_result.returncode != 0:
		raise RuntimeError("Liftoff did not run successfully.")

	subprocess.run(["rm", str(folder / "reference.fa.mmi")])

	gff3_with_version = Gff3Parser.read_gff3_as_bytes_add_isoformcheck_version_to_start(tmp_folder / "reference_annotation.gff3", filtered_refannotation_hash)
	with open(str(sample_annotation_folder / "reference.gff3.gz"), "wb") as compressed_gff3:
		compress_command = ["gzip"]
		gzip_result = subprocess.run(compress_command, input=gff3_with_version, stdout=compressed_gff3)
		if gzip_result.returncode != 0:
			raise RuntimeError("gzip did not run successfully")

	agc_command = [agc_path, "create", str(folder / "reference.fa")]
	with open(str(folder / "sequences.agc"), "wb") as agc_file:
		agc_result = subprocess.run(agc_command, stdout=agc_file)
		if agc_result.returncode != 0:
			raise RuntimeError("agc did not run successfully")

	initialize_sqlite_db_schema(folder / "sample_info.db")
	insert_reference_annotations_to_db(folder / "sample_info.db", folder / "reference.gff3")

	shutil.rmtree(tmp_folder)

	HandleAssembly.add_sample_proteins_to_database(folder, sample_annotation_folder / "reference.gff3.gz", ref_fasta, "reference", "reference")
