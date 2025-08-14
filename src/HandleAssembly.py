#!/usr/bin/env python

import os
import sys
import gzip
import pathlib
import subprocess
import shutil
import sqlite3
import datetime
import Gff3Parser
import TranscriptExtractor
import SequenceReader

def prepare_fasta(input_path, output_path):
	"""
	Creates an unzipped fasta out of a sequence file. Input may be fasta or fastq, possibly gzipped.

	Args:
		input_path: Path to input sequence in fasta/fastq/gzipped format
		output_path: Path to resulting unzipped fasta
	"""

	with open(output_path, "w") as f:
		for (name, sequence) in SequenceReader.stream_sequences(input_path):
			print(">" + name, file=f)
			print(sequence, file=f)

def add_sample_proteins_to_database(database_file, sample_annotation, sample_fasta, sample_name, sample_haplotype):
	"""
	Adds proteins of a sample to the database. Assumes that liftoff has already been ran and sample annotation is in subfolder sample_annotation

	Args:
		database_file: Base folder
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
	"""

	print(f"step 1 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	sample_transcripts = TranscriptExtractor.process_sample_transcripts(sample_fasta, sample_annotation)

	print(f"step 2 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	sample_contig_lens = {}
	for name, sequence in SequenceReader.stream_sequences(sample_fasta):
		sample_contig_lens[name] = len(sequence)

	print(f"step 3 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	(gene_locations, transcript_locations) = Gff3Parser.parse_gene_transcript_locations(sample_annotation)

	print(f"step 4 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	processed_transcripts = []
	for transcript_id, sequence, extra_copy in sample_transcripts:
		name = transcript_id
		if extra_copy != 0:
			name = "_".join(name.split("_")[:-1])
		processed_transcripts.append((name, sequence, transcript_locations[transcript_id]))

	transcript_id_map = {}
	allele_name_map = {}
	sample_contig_db_ids = {}
	novel_alleles = set()
	with sqlite3.connect(str(database_file)) as connection:
		cursor = connection.cursor()
		cursor.arraysize = 10000
		print(f"step 5 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
		if sample_id is None:
			cursor.execute("INSERT INTO Sample (Name) VALUES (?)", (sample_name,))
			sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
			assert sample_id is not None
		sample_id = sample_id[0]
		print(f"step 6 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		cursor.execute("INSERT INTO Haplotype (SampleId, Haplotype) VALUES (?, ?)", (sample_id, sample_haplotype))
		print(f"step 7 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		haplotype_id = cursor.execute("SELECT Id FROM Haplotype WHERE SampleId=? AND Haplotype=?", (sample_id, sample_haplotype)).fetchone()
		assert haplotype_id is not None
		haplotype_id = haplotype_id[0]
		print(f"step 8 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for name, length in sample_contig_lens.items():
			cursor.execute("INSERT INTO SampleContig (HaplotypeId, Name, Length) VALUES (?, ?, ?)", (haplotype_id, name, length))
		print(f"step 9 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for db_id, name in cursor.execute("SELECT Id, Name FROM SampleContig WHERE HaplotypeId=?", (haplotype_id,)):
			sample_contig_db_ids[name] = db_id
		print(f"step 10 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for db_id, transcript_id in cursor.execute("SELECT Id, Name FROM Transcript").fetchall():
			transcript_id_map[transcript_id] = db_id
		print(f"step 11 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for allele_id, transcript_id, transcript_sequence in cursor.execute("SELECT Id, TranscriptId, Sequence FROM Allele").fetchall():
			allele_name_map[(transcript_id, transcript_sequence)] = allele_id
		print(f"step 12 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for (transcript_id, sequence, location) in processed_transcripts:
			transcript_db_id = transcript_id_map[transcript_id]
			if (transcript_db_id, sequence) not in allele_name_map: novel_alleles.add((transcript_db_id, sequence))
		print(f"step 13 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		if len(novel_alleles) >= 1:
			print(f"step 14 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			cursor.executemany("INSERT INTO Allele (TranscriptId, Sequence) VALUES (?, ?)", novel_alleles)
		print(f"step 14.5 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		if len(novel_alleles) >= 1:
			for allele_id, transcript_id, transcript_sequence in cursor.execute("SELECT Allele.Id, Transcript.Id, Allele.Sequence FROM Transcript INNER JOIN Allele ON Allele.TranscriptId = Transcript.Id").fetchall():
				allele_name_map[(transcript_id, transcript_sequence)] = allele_id
		print(f"step 15 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		insert_lines = []
		for transcript_id, sequence, location in processed_transcripts:
			insert_lines.append((haplotype_id, allele_name_map[(transcript_id_map[transcript_id], sequence)], sample_contig_db_ids[location[0]], location[1], location[2], location[3]))
		print(f"step 16 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		cursor.executemany("INSERT INTO SampleProtein (HaplotypeId, AlleleId, SampleContigId, SampleLocationStrand, SampleLocationStart, SampleLocationEnd) VALUES (?, ?, ?, ?, ?, ?)", insert_lines)
		print(f"step 17 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		connection.commit()
		print(f"step 18 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	print(f"step 19 time {datetime.datetime.now().astimezone()}", file=sys.stderr)

def handle_new_sample_liftoff_and_transcripts(database_folder, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path, agc_path):
	"""
	Handles a new sample. Adds the sample sequences to the agc database, runs liftoff and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_sequence: Sequence file of sample in fasta/fastq/gzip format
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
	"""
	sample_exists = check_if_sample_exists(database_folder, sample_name, sample_haplotype)
	if sample_exists:
		raise RuntimeError(f"Sample already exists: name \"{sample_name}\" haplotype \"{sample_haplotype}\"")

	tmp_folder = database_folder / ("tmp_" + sample_name + "_" + sample_haplotype)
	os.makedirs(tmp_folder, exist_ok=False)
	try:
		handle_new_sample_liftoff_use_tmp_folder(database_folder, tmp_folder, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path, agc_path)
		add_sample_proteins_to_database(database_folder / "sample_info.db", database_folder / "sample_annotations" / (sample_name + "_" + sample_haplotype + ".gff3.gz"), tmp_folder / (sample_name + "_" + sample_haplotype + ".fa"), sample_name, sample_haplotype)
	finally:
		shutil.rmtree(tmp_folder)

def check_if_sample_exists(database_folder, sample_name, sample_haplotype):
	"""
	Check if a sample haplotype exists.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
	"""
	with sqlite3.connect(str(database_folder / "sample_info.db")) as connection:
		cursor = connection.cursor()
		sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
		if sample_id is None: return False
		sample_id = sample_id[0]
		haplotype_id = cursor.execute("SELECT Id FROM Haplotype WHERE SampleId=? AND Haplotype=?", (sample_id, sample_haplotype)).fetchone()
		if haplotype_id is not None: return True
	return False

def handle_new_sample_liftoff_use_tmp_folder(database_folder, tmp_folder, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path, agc_path):
	"""
	Handles a new sample liftoff. Adds the sample sequences to the agc database, runs liftoff and stores result in sample_annotation folder.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		tmp_folder: Folder where temporary files will be stored.
		sample_sequence: Sequence file of sample in fasta/fastq/gzip format
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
	"""

	sample_annotation_folder = database_folder / "sample_annotations"
	reference_sequence_path = database_folder / "reference.fa"
	temp_file_path = tmp_folder / (sample_name + "_" + sample_haplotype + ".fa")
	prepare_fasta(sample_sequence, str(temp_file_path))

	agc_command = [agc_path, "append", str(database_folder / "sequences.agc"), str(temp_file_path)]
	with open(str(tmp_folder / "tmp.agc"), "wb") as new_agc:
		agc_result = subprocess.run(agc_command, stdout=new_agc)
		if agc_result.returncode != 0:
			raise RuntimeError("agc did not run successfully.")

	liftoff_command = [liftoff_path, "-db", str(database_folder / "reference.gff3_db"), "-o", str(tmp_folder / "tmp_annotation.gff3"), "-p", str(num_threads), "-sc", "0.95", "-copies", "-polish", "-dir", str(tmp_folder / "intermediate_files"), "-u", str(tmp_folder / "unmapped_features.txt"), str(temp_file_path), str(reference_sequence_path)]
	print(f"Running liftoff with command:", file=sys.stderr)
	print(f"{' '.join(liftoff_command)}", file=sys.stderr)
	liftoff_result = subprocess.run(liftoff_command)
	if liftoff_result.returncode != 0:
		raise RuntimeError("Liftoff did not run successfully.")

	gzip_command = ["gzip"]
	with open(tmp_folder / "tmp_annotation.gff3") as raw_gff3:
		with open(sample_annotation_folder / (sample_name + "_" + sample_haplotype + ".gff3.gz"), "wb") as compressed_gff3:
			gzip_result = subprocess.run(gzip_command, stdin=raw_gff3, stdout=compressed_gff3)
			if gzip_result.returncode != 0:
				raise RuntimeError("gzip did not run successfully")

	move_command = ["mv", str(tmp_folder / "tmp.agc"), str(database_folder / "sequences.agc")]
	move_result = subprocess.run(move_command)
	if move_result.returncode != 0:
		raise RuntimeError("mv did not run successfully")

def check_if_sample_sex_chromosome_annotations_are_fine(database_file, sample_name):
	"""
	Checks if sample's gene annotations in sex chromosomes are as expected. Sample should not have both chrX and chrY genes in the sample haplotype, and sample should not have chrY genes in mat haplotype.

	Args:
		database_file: Path to sample info sql db
		sample_name: Name of sample

	Returns:
		List with strings describing what's wrong. If all OK then list is empty
	"""

	select_command = \
	"""
		SELECT DISTINCT Gene.ReferenceChromosome, Haplotype.Haplotype
		FROM SampleProtein
		INNER JOIN Allele ON Allele.Id = SampleProtein.AlleleId
		INNER JOIN Transcript ON Transcript.Id = Allele.TranscriptId
		INNER JOIN Gene ON Gene.Id = Transcript.GeneId
		INNER JOIN Haplotype ON Haplotype.Id = SampleProtein.HaplotypeId
		INNER JOIN Sample ON Haplotype.SampleId = Sample.Id
		WHERE Sample.Name=?
	"""
	result = []
	with sqlite3.connect(str(database_file)) as connection:
		cursor = connection.cursor()
		reference_chromosomes_per_haplotype = {}
		for reference_chromosome, haplotype in cursor.execute(select_command, (sample_name,)).fetchall():
			if haplotype not in reference_chromosomes_per_haplotype: reference_chromosomes_per_haplotype[haplotype] = set()
			reference_chromosomes_per_haplotype[haplotype].add(reference_chromosome)
		for haplotype in reference_chromosomes_per_haplotype:
			if haplotype == "mat" and "chrY" in reference_chromosomes_per_haplotype[haplotype]:
				result.append("Haplotype 'mat' contains genes from chromosome Y")
			if "chrX" in reference_chromosomes_per_haplotype[haplotype] and "chrY" in reference_chromosomes_per_haplotype[haplotype]:
				result.append(f"Haplotype '{haplotype}' contains genes from both chromosome X and chromosome Y")
	return result
