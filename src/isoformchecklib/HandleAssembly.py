#!/usr/bin/env python

import os
import sys
import gzip
import pathlib
import subprocess
import shutil
import sqlite3
import datetime
import threading
import queue
import isoformchecklib.Gff3Parser as Gff3Parser
import isoformchecklib.TranscriptExtractor as TranscriptExtractor
import isoformchecklib.SequenceReader as SequenceReader
import isoformchecklib.Util as Util
import isoformchecklib.DatabaseOperations as DatabaseOperations

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

def validate_gff3_is_same_isoformcheck_version(annotation_file, refannotation_hash):
	"""
	Validates that a gff3 file has the same IsoformCheck DB version and reference annotation hash as currently.

	Args:
		annotation_file: Path to input gff3
		refannotation_hash: String with reference annotation hash

	Returns:
		True if both IsoformCheck DB version and reference annotation hash match exactly. False otherwise
	"""
	wanted_version_string = "# IsoformCheck DB version " + Util.DBVersion + " annotation hash " + refannotation_hash
	with Util.open_maybe_gzipped(annotation_file) as f:
		for l in f:
			if l.strip() == wanted_version_string:
				return True
			if not l.startswith('#'):
				return False
	return False

def compare_samples_to_database(base_folder, sample_table_file, num_threads, liftoff_path, force):
	"""
	Compares multiple new samples to the database. Runs liftoff if necessary and gets transcripts, finds novel isoforms and allele sets, and sample allele sets. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		base_folder: Folder where database files are located. Type should be pathlib.Path
		sample_table_file: Tsv file with sample info
		num_threads: Number of threads
		liftoff_path: Path to liftoff executable
		force: Force add samples even if annotation version does not match

	Returns:
		tuple of (novel_isoforms, novel_allelesets, sample_allelesets)
		Format of novel_isoforms is [(transcript, isoform_name, isoform_sequence)]
		Format of novel_allelesets is [(transcript, sample, alleleset)]
		Format of sample_allelesets is [(transcript, sample, alleleset)]
	"""
	refannotation_hash = Util.get_refannotation_hash(base_folder / "info.txt")
	sample_info = read_sample_info_from_table(sample_table_file)
	dupes = check_sample_duplicates(sample_info)
	if dupes:
		raise Util.ParameterError(f"Duplicate sample: sample name \"{dupes[0]}\" haplotype \"{dupes[1]}\"")
	incongruents = check_samples_with_incongruent_haplotypes(sample_info)
	if len(incongruents) > 0:
		invalid_samples = []
		for sample, haplotypes in incongruents:
			invalid_samples += f"sample \"{sample}\" haplotypes \"{",".join(haplotypes)}\""
		invalid_samples_string = " ; ".join(invalid_samples)
		raise Util.ParameterError(f"Invalid haplotypes in samples. All samples should have two haplotypes, either \"1\" and \"2\" or \"mat\" and \"pat\". Invalid samples and their haplotypes: {invalid_samples_string}")
	for sample_name, sample_haplotype, _, _ in sample_info:
		sample_exists = check_if_sample_exists(base_folder, sample_name, sample_haplotype)
		if sample_exists:
			raise Util.ParameterError(f"Sample already exists in the database: name \"{sample_name}\" haplotype \"{sample_haplotype}\"")
	tmp_prefix = Util.make_random_prefix()
	tmp_base_folder = base_folder / ("tmp_" + tmp_prefix)
	result = None
	try:
		os.makedirs(tmp_base_folder, exist_ok=False)
		annotation_folder = tmp_base_folder / "sample_annotations"
		os.makedirs(annotation_folder, exist_ok=False)
		sample_info_with_annotations = []
		for sample_name, sample_haplotype, sample_sequence, annotation in sample_info:
			if annotation:
				Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Using annotations for sample \"{sample_name}\" haplotype \"{sample_haplotype}\" from path \"{annotation}\"", file=sys.stderr)
				sample_info_with_annotations.append((sample_name, sample_haplotype, sample_sequence, annotation))
			else:
				Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Running liftoff for sample \"{sample_name}\" haplotype \"{sample_haplotype}\"", file=sys.stderr)
				sample_annotation_file = annotation_folder / (sample_name + "_" + sample_haplotype + ".gff3.gz")
				tmp_folder = tmp_base_folder / ("tmp_" + sample_name + "_" + sample_haplotype)
				os.makedirs(tmp_folder, exist_ok=False)
				handle_new_sample_liftoff_use_tmp_folder(base_folder, tmp_folder, sample_annotation_file, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path)
				sample_info_with_annotations.append((sample_name, sample_haplotype, tmp_folder / (sample_name + "_" + sample_haplotype + ".fa"), sample_annotation_file))
		samples_with_mismatch = []
		for name, hap, fasta_file, annotation_file in sample_info_with_annotations:
			version_match = validate_gff3_is_same_isoformcheck_version(annotation_file, refannotation_hash)
			if not version_match:
				samples_with_mismatch.append((name, hap))
		if len(samples_with_mismatch) >= 1:
			if not force:
				raise Util.ParameterError("Sample annotation does not match IsoformCheck version. Rerun liftover for samples, or if you are sure about what you are doing you can force insert with --force. Samples and haplotypes with invalid versions: " + ", ".join(name + " " + hap for name, hap in samples_with_mismatch))
			else:
				print(f"{datetime.datetime.now().astimezone()}: Sample annotation does not match IsoformCheck version. Adding samples anyway due to --force. Samples and haplotypes with invalid versions: " + ", ".join(name + " " + hap for name, hap in samples_with_mismatch), file=sys.stderr)
		result = get_novel_isoforms_allelesets(base_folder, sample_info_with_annotations, num_threads)
	finally:
		shutil.rmtree(tmp_base_folder)
	return result

def get_sample_transcripts_and_contig_lens(sample_info, num_threads):
	"""
	Reads sample annotations and sequences, processes transcripts and contig lens. Assumes annotations are already present.

	Args:
		sample_info: List of [(sample_name, sample_haplotype, sample_sequence_file_path, sample_annotation_file_path)]
		num_threads: Number of threads

	Returns:
		(sample_transcripts, sample_contig_lengths)
		Format of sample_transcripts is [(transcript_id, sequence, transcript_location)]
		Format of sample_contig_lengths is dict of (contig -> length)
	"""
	threads = []
	def process_transcripts(sample_info, input_id_queue, output_transcript_queue):
		while True:
			if input_id_queue.empty(): return
			index = input_id_queue.get()
			sample_name, sample_haplotype, sample_fasta, sample_annotation = sample_info[index]
			Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Get transcripts of sample {sample_name} haplotype {sample_haplotype}", file=sys.stderr)
			result_here = []
			(sample_transcripts, gene_locations, contig_lengths) = TranscriptExtractor.process_sample_transcripts_and_contigs(sample_fasta, sample_annotation)
			for transcript_id, sequence, extra_copy, location in sample_transcripts:
				result_here.append((transcript_id, sequence, location))
			output_transcript_queue.put((index, result_here, contig_lengths))
			Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Got transcripts of sample {sample_name} haplotype {sample_haplotype}", file=sys.stderr)
	input_ids = queue.Queue(len(sample_info))
	output_transcripts = queue.Queue(len(sample_info))
	for i in range(0, len(sample_info)):
		input_ids.put(i)
	for i in range(0, num_threads):
		threads.append(threading.Thread(target=process_transcripts, args=(sample_info, input_ids, output_transcripts)))
	for i in range(0, num_threads):
		threads[i].start()
	for i in range(0, num_threads):
		threads[i].join()
	Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Merging sample transcripts and contig lengths", file=sys.stderr)
	all_processed_transcripts = []
	all_sample_contig_lens = []
	for i in range(0, len(sample_info)):
		all_processed_transcripts.append(None)
		all_sample_contig_lens.append(None)
	for i in range(0, len(sample_info)):
		index, processed_transcripts, sample_contig_lens = output_transcripts.get()
		assert index < len(all_processed_transcripts)
		assert all_processed_transcripts[index] is None
		all_processed_transcripts[index] = processed_transcripts
		assert all_sample_contig_lens[index] is None
		all_sample_contig_lens[index] = sample_contig_lens
	for i in range(0, len(sample_info)):
		assert all_processed_transcripts[i] is not None
		assert all_sample_contig_lens[i] is not None
	return all_processed_transcripts, all_sample_contig_lens

def get_novel_isoforms_allelesets(base_folder, sample_info, num_threads):
	"""
	Gets the novel isoforms of samples, novel allelesets, and all sample allelesets. Assumes that liftoff has already been ran

	Args:
		base_folder: Base folder
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_fasta_file_path, sample_annotation_file_path)
		num_threads: Number of threads to use. Processing each sample uses one thread but multiple samples can be processed in parallel

	Returns:
		tuple of (novel_isoforms, novel_allelesets, sample_allelesets)
		Format of novel_isoforms is [(transcript, isoform_name, isoform_sequence)]
		Format of novel_allelesets is [(transcript, sample, alleleset)]
		Format of sample_allelesets is [(transcript, sample, alleleset)]
	"""
	Util.verbose_print(2, f"step 1 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	all_processed_transcripts, all_sample_contig_lens = get_sample_transcripts_and_contig_lens(sample_info, num_threads)
	Util.verbose_print(2, f"step 2 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	transcript_id_map = {}
	isoform_name_map = {}
	isoform_sequences = DatabaseOperations.get_isoform_sequences_from_file(str(base_folder / "isoforms.fa"))
	novel_isoforms = set()
	next_novel_per_transcript = {}
	novel_sample_transcript_alleles = {}
	all_transcripts = set()
	with sqlite3.connect(str(base_folder / "sample_info.db")) as connection:
		cursor = connection.cursor()
		cursor.arraysize = 10000
		for db_id, transcript_id in cursor.execute("SELECT Id, Name FROM Transcript").fetchall():
			transcript_id_map[transcript_id] = db_id
		for isoform_name, transcript_id, transcript_sequence_id in cursor.execute("SELECT Name, TranscriptId, SequenceId FROM Isoform").fetchall():
			isoform_name_map[(transcript_id, isoform_sequences[transcript_sequence_id])] = isoform_name
		for sampleindex in range(0, len(sample_info)):
			sample_name = sample_info[sampleindex][0]
			processed_transcripts = all_processed_transcripts[sampleindex]
			for (transcript_id, sequence, location) in processed_transcripts:
				transcript_db_id = transcript_id_map[transcript_id]
				if (transcript_db_id, sequence) not in isoform_name_map:
					if transcript_db_id not in next_novel_per_transcript: next_novel_per_transcript[transcript_db_id] = 0
					novel_name = "novel_" + DatabaseOperations.get_isoform_name(next_novel_per_transcript[transcript_db_id])
					novel_isoforms.add((transcript_id, novel_name, sequence))
					isoform_name_map[transcript_db_id, sequence] = novel_name
					next_novel_per_transcript[transcript_db_id] += 1
				isoform_name = isoform_name_map[transcript_db_id, sequence]
				if sample_name not in novel_sample_transcript_alleles: novel_sample_transcript_alleles[sample_name] = {}
				if transcript_id not in novel_sample_transcript_alleles[sample_name]: novel_sample_transcript_alleles[sample_name][transcript_id] = []
				novel_sample_transcript_alleles[sample_name][transcript_id].append(isoform_name)
	for sample in novel_sample_transcript_alleles:
		for transcript in transcript_id_map:
			alleleset = ()
			if transcript in novel_sample_transcript_alleles[sample]:
				novel_sample_transcript_alleles[sample][transcript] = list(novel_sample_transcript_alleles[sample][transcript])
				novel_sample_transcript_alleles[sample][transcript].sort(key=lambda x: DatabaseOperations.isoform_sort_order(x))
				alleleset = tuple(novel_sample_transcript_alleles[sample][transcript])
			novel_sample_transcript_alleles[sample][transcript] = alleleset
		for transcript_id in transcript_id_map:
			if transcript_id in novel_sample_transcript_alleles[sample]: continue
			alleleset = ()
			novel_sample_transcript_alleles[sample][transcript] = alleleset
	existing_allelesets = {}
	with sqlite3.connect(str(base_folder / "sample_info.db")) as connection:
		cursor = connection.cursor()
		cursor.arraysize = 10000
		sample_transcript_alleles = {}
		for isoform_name, transcript, sample_name in cursor.execute("SELECT Isoform.Name, Transcript.Name, Sample.Name FROM SampleAllele INNER JOIN Isoform ON SampleAllele.IsoformId = Isoform.Id INNER JOIN Transcript ON Isoform.TranscriptId = Transcript.Id INNER JOIN Haplotype ON SampleAllele.HaplotypeId = Haplotype.Id INNER JOIN Sample ON Sample.Id = Haplotype.SampleId"):
			if sample_name not in sample_transcript_alleles: sample_transcript_alleles[sample_name] = {}
			if transcript not in sample_transcript_alleles[sample_name]: sample_transcript_alleles[sample_name][transcript] = []
			sample_transcript_alleles[sample_name][transcript].append(isoform_name)
		for sample in sample_transcript_alleles:
			for transcript in transcript_id_map:
				alleleset = ()
				if transcript in sample_transcript_alleles[sample]:
					sample_transcript_alleles[sample][transcript].sort(key=lambda x: DatabaseOperations.isoform_sort_order(x))
					alleleset = tuple(sample_transcript_alleles[sample][transcript])
				if transcript not in existing_allelesets: existing_allelesets[transcript] = set()
				existing_allelesets[transcript].add(alleleset)
	novel_allelesets = []
	all_sample_allelesets = []
	for sample in novel_sample_transcript_alleles:
		for transcript in novel_sample_transcript_alleles[sample]:
			if novel_sample_transcript_alleles[sample][transcript] not in existing_allelesets[transcript]:
				novel_allelesets.append((transcript, sample, novel_sample_transcript_alleles[sample][transcript]))
			all_sample_allelesets.append((transcript, sample, novel_sample_transcript_alleles[sample][transcript]))
	novel_isoforms = list(novel_isoforms)
	novel_isoforms.sort(key=lambda x: (x[0], x[1], DatabaseOperations.isoform_sort_order(x[2])))
	novel_allelesets.sort(key=lambda x: (x[0], x[1], DatabaseOperations.alleleset_sort_order(x[2])))
	all_sample_allelesets.sort(key=lambda x: (x[0], x[1], DatabaseOperations.alleleset_sort_order(x[2])))
	return (novel_isoforms, novel_allelesets, all_sample_allelesets)

def add_multiple_sample_proteins_to_database(base_folder, sample_info, num_threads):
	"""
	Adds proteins of a sample to the database. Assumes that liftoff has already been ran

	Args:
		base_folder: Base folder
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_fasta_file_path, sample_annotation_file_path)
		num_threads: Number of threads to use. Processing each sample uses one thread but multiple samples can be processed in parallel
	"""
	Util.verbose_print(2, f"step 1 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	all_processed_transcripts, all_sample_contig_lens = get_sample_transcripts_and_contig_lens(sample_info, num_threads)
	Util.verbose_print(2, f"step 2 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	transcript_id_map = {}
	isoform_name_map = {}
	sample_contig_db_ids = []
	all_haplotype_ids = []
	isoform_sequences = DatabaseOperations.get_isoform_sequences_from_file(str(base_folder / "isoforms.fa"))
	novel_isoforms = set()
	with sqlite3.connect(str(base_folder / "sample_info.db")) as connection:
		cursor = connection.cursor()
		cursor.arraysize = 10000
		Util.verbose_print(2, f"step 3 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for db_id, transcript_id in cursor.execute("SELECT Id, Name FROM Transcript").fetchall():
			transcript_id_map[transcript_id] = db_id
		Util.verbose_print(2, f"step 4 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for isoform_id, transcript_id, transcript_sequence_id in cursor.execute("SELECT Id, TranscriptId, SequenceId FROM Isoform").fetchall():
			isoform_name_map[(transcript_id, isoform_sequences[transcript_sequence_id])] = isoform_id
		Util.verbose_print(2, f"step 5 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for sampleindex in range(0, len(sample_info)):
			sample_name = sample_info[sampleindex][0]
			sample_haplotype = sample_info[sampleindex][1]
			Util.verbose_print(2, f"step 6 sample {sample_name} hap {sample_haplotype} time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			sample_contig_db_ids.append({})
			sample_contig_lens = all_sample_contig_lens[sampleindex]
			processed_transcripts = all_processed_transcripts[sampleindex]
			sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
			if sample_id is None:
				cursor.execute("INSERT INTO Sample (Name) VALUES (?)", (sample_name,))
				sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
				assert sample_id is not None
			sample_id = sample_id[0]
			Util.verbose_print(2, f"step 7 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			cursor.execute("INSERT INTO Haplotype (SampleId, Haplotype) VALUES (?, ?)", (sample_id, sample_haplotype))
			Util.verbose_print(2, f"step 8 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			haplotype_id = cursor.execute("SELECT Id FROM Haplotype WHERE SampleId=? AND Haplotype=?", (sample_id, sample_haplotype)).fetchone()
			assert haplotype_id is not None
			haplotype_id = haplotype_id[0]
			all_haplotype_ids.append(haplotype_id)
			Util.verbose_print(2, f"step 9 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for name, length in sample_contig_lens.items():
				cursor.execute("INSERT INTO SampleContig (HaplotypeId, Name, Length) VALUES (?, ?, ?)", (haplotype_id, name, length))
			Util.verbose_print(2, f"step 10 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for db_id, name in cursor.execute("SELECT Id, Name FROM SampleContig WHERE HaplotypeId=?", (haplotype_id,)):
				sample_contig_db_ids[sampleindex][name] = db_id
			Util.verbose_print(2, f"step 11 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for (transcript_id, sequence, location) in processed_transcripts:
				transcript_db_id = transcript_id_map[transcript_id]
				if (transcript_db_id, sequence) not in isoform_name_map: novel_isoforms.add((transcript_db_id, sequence))
		Util.verbose_print(2, f"step 12 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		if len(novel_isoforms) >= 1:
			inserted_isoforms = DatabaseOperations.add_isoforms_to_fasta(base_folder / "isoforms.fa", novel_isoforms)
			Util.verbose_print(2, f"step 13 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			cursor.executemany("INSERT INTO isoform (TranscriptId, SequenceId) VALUES (?, ?)", inserted_isoforms)
			Util.verbose_print(2, f"step 13.5 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			isoform_sequences = DatabaseOperations.get_isoform_sequences_from_file(str(base_folder / "isoforms.fa"))
			Util.verbose_print(2, f"step 14 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for isoform_id, transcript_id, transcript_sequence_id in cursor.execute("SELECT isoform.Id, Transcript.Id, isoform.SequenceId FROM Transcript INNER JOIN isoform ON isoform.TranscriptId = Transcript.Id").fetchall():
				isoform_name_map[(transcript_id, isoform_sequences[transcript_sequence_id])] = isoform_id
		Util.verbose_print(2, f"step 15 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		insert_lines = []
		for sampleindex in range(0, len(sample_info)):
			processed_transcripts = all_processed_transcripts[sampleindex]
			for transcript_id, sequence, location in processed_transcripts:
				insert_lines.append((all_haplotype_ids[sampleindex], isoform_name_map[(transcript_id_map[transcript_id], sequence)], sample_contig_db_ids[sampleindex][location[0]], location[1], location[2], location[3]))
		Util.verbose_print(2, f"step 16 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		cursor.executemany("INSERT INTO SampleAllele (HaplotypeId, IsoformId, SampleContigId, SampleLocationStrand, SampleLocationStart, SampleLocationEnd) VALUES (?, ?, ?, ?, ?, ?)", insert_lines)
		Util.verbose_print(2, f"step 17 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		connection.commit()
		Util.verbose_print(2, f"step 18 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	Util.verbose_print(2, f"step 19 time {datetime.datetime.now().astimezone()}", file=sys.stderr)

def add_sample_proteins_to_database(base_folder, sample_annotation, sample_fasta, sample_name, sample_haplotype):
	"""
	Adds proteins of a sample to the database. Assumes that liftoff has already been ran

	Args:
		base_folder: Base folder
		sample_annotation: Path of sample gff3 annotation file
		sample_fasta: Path of sample fasta file
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
	"""
	add_multiple_sample_proteins_to_database(base_folder, [(sample_name, sample_haplotype, sample_fasta, sample_annotation)], 1)

def handle_one_new_sample_liftoff_and_transcripts(base_folder, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path, agc_path):
	"""
	Handles one new sample. Adds the sample sequences to the agc database, runs liftoff and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		base_folder: Folder where database files are located. Type should be pathlib.Path
		sample_sequence: Sequence file of sample in fasta/fastq/gzip format
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
	"""
	sample_info = [(sample_name, sample_haplotype, sample_sequence, None)]
	handle_multiple_new_samples_liftoff_and_transcripts(base_folder, sample_info, num_threads, liftoff_path, agc_path, False)

def run_liftoff(database_folder, input_sequence, output, num_threads, liftoff_path):
	"""
	Runs liftoff for one haplotype. Input should be a file containing the contigs / scaffolds of a single haplotype.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		input_sequence: Path to input sequence file
		output: Path to output annotation file. Should have ending .gff3.gz
		liftoff_path: Path to liftoff executable
	"""
	sample_name = Util.make_random_prefix()
	sample_haplotype = Util.make_random_prefix()
	tmp_folder = database_folder / ("tmp_" + sample_name + "_" + sample_haplotype)
	try:
		os.makedirs(tmp_folder, exist_ok=False)
		handle_new_sample_liftoff_use_tmp_folder(database_folder, tmp_folder, output, input_sequence, sample_name, sample_haplotype, num_threads, liftoff_path)
	finally:
		shutil.rmtree(tmp_folder)

def read_sample_info_from_table(sample_table_file):
	"""
	Reads sample info table from tsv file

	Args:
		sample_table_file: Tsv file with descriptions of the samples

	Returns:
		[(sample_name, sample_haplotype, sample_assembly_file_path, sample_annotation_file_path)]
	"""
	sample_info = []
	row_number = 0
	with open(sample_table_file) as f:
		has_annotation = None
		for l in f:
			if len(l.strip()) == 0: continue
			row_number += 1
			parts = l.strip().split("\t")
			if len(parts) != 3 and len(parts) != 4:
				raise Util.ParameterError(f"Sample table file has wrong format. File should be a tab-separated file with three or four columns per row. Row {row_number} has {len(parts)} columns.")
			if row_number == 1:
				if parts[0].lower() != "sample":
					raise Util.ParameterError(f"Sample table file has wrong format. First column should be sample name. Header instead has \"{parts[0]}\".")
				if parts[1].lower() != "haplotype":
					raise Util.ParameterError(f"Sample table file has wrong format. Second column should be haplotype name. Header instead has \"{parts[1]}\".")
				if parts[2].lower() != "assembly":
					raise Util.ParameterError(f"Sample table file has wrong format. Third column should be assembly file. Header instead has \"{parts[2]}\".")
				if len(parts) >= 4 and parts[3].lower() != "annotation":
					raise Util.ParameterError(f"Sample table file has wrong format. Fourth column should be annotation file. Header instead has \"{parts[3]}\".")
				has_annotation = (len(parts) == 4)
				continue
			if parts[1] not in ["1", "2", "mat", "pat"]:
				raise Util.ParameterError(f"Sample table file has wrong format. Row {row_number} haplotype is \"{parts[2]}\", expected one of \"1\", \"2\", \"mat\", \"pat\"")
			annotation = None
			if has_annotation:
				if len(parts) >= 4:
					if parts[3] is not None and len(parts[3]) >= 1:
						annotation = parts[3]
			sample_info.append((parts[0], parts[1], parts[2], annotation))
			if not Util.file_exists(parts[2]):
				raise Util.ParameterError(f"Sample assembly file in row {row_number} cannot be read: file \"{parts[2]}\", sample \"{parts[0]}\" haplotype \"{parts[1]}\"")
			if annotation is not None and not Util.file_exists(annotation):
				raise Util.ParameterError(f"Sample annotation file in row {row_number} cannot be read: file \"{parts[3]}\", sample \"{parts[0]}\" haplotype \"{parts[1]}\"")
	return sample_info

def check_sample_duplicates(sample_info):
	"""
	Checks if sample info has duplicates

	Args:
		sample_info: [(sample_name, sample_haplotype, sample_assembly_file_path, sample_annotation_file_path)]

	Returns:
		None if no duplicates are present, (sample, haplotype) of one duplicate otherwise
	"""
	has_duplicates = False
	sample_name_and_haplotypes = set()
	for row in sample_info:
		key = (row[0], row[1])
		if key in sample_name_and_haplotypes:
			return key
		sample_name_and_haplotypes.add(key)

def check_samples_with_incongruent_haplotypes(sample_info):
	"""
	Checks if every sample has two haplotypes, which are named either 1 and 2 or mat and pat

	Args:
		sample_info: [(sample_name, sample_haplotype, sample_assembly_file_path, sample_annotation_file_path)]

	Returns:
		List of samples with wrong haplotypes. Empty list if all fine. Format is [(sample_name, [sample_haplotype])]
	"""
	sample_haplotypes = {}
	for row in sample_info:
		sample_name = row[0]
		sample_haplotype = row[1]
		if sample_name not in sample_haplotypes: sample_haplotypes[sample_name] = []
		sample_haplotypes[sample_name].append(sample_haplotype)
	incongruent = []
	for sample in sample_haplotypes:
		sample_haplotypes[sample].sort()
		if len(sample_haplotypes[sample]) == 2:
			if sample_haplotypes[sample][0] == "1" and sample_haplotypes[sample][1] == "2":
				continue
			if sample_haplotypes[sample][0] == "mat" and sample_haplotypes[sample][1] == "pat":
				continue
		incongruent.append((sample, sample_haplotypes[sample]))
	return incongruent

def handle_multiple_new_samples_liftoff_and_transcripts_from_table(database_folder, sample_table_file, num_threads, liftoff_path, agc_path, force):
	"""
	Handles multiple samples. Adds the sample sequences to the agc database, runs liftoff and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.
	If sample table file has a column for annotation, then annotations are copied from there instead of rerun

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_table_file: Tsv file with descriptions of the samples
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
		force: Force add samples even if annotation version does not match
	"""
	sample_info = read_sample_info_from_table(sample_table_file)
	has_dupes = check_sample_duplicates(sample_info)
	if has_dupes:
		raise Util.ParameterError(f"Duplicate sample and haplotype: sample \"{has_dupes[0]}\" haplotype \"{has_dupes[1]}\"")
	Util.verbose_print(0, f"{datetime.datetime.now().astimezone()}: Adding {len(sample_info)} new assemblies", file=sys.stderr)
	handle_multiple_new_samples_liftoff_and_transcripts(database_folder, sample_info, num_threads, liftoff_path, agc_path, force)
	Util.verbose_print(0, f"{datetime.datetime.now().astimezone()}: Finished adding {len(sample_info)} new assemblies", file=sys.stderr)

def handle_multiple_new_samples_liftoff_and_transcripts(base_folder, sample_info, num_threads, liftoff_path, agc_path, force):
	"""
	Handles multiple new samples. Adds the sample sequences to the agc database, runs liftoff if necessary and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		base_folder: Folder where database files are located. Type should be pathlib.Path
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_sequence_file_path, sample_annotation_file_path)
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
		force: Force add samples even if annotation version does not match
	"""
	refannotation_hash = Util.get_refannotation_hash(base_folder / "info.txt")
	for sample_name, sample_haplotype, _, _ in sample_info:
		sample_exists = check_if_sample_exists(base_folder, sample_name, sample_haplotype)
		if sample_exists:
			raise Util.ParameterError(f"Sample already exists: name \"{sample_name}\" haplotype \"{sample_haplotype}\"")

	tmp_prefix = Util.make_random_prefix()
	tmp_base_folder = base_folder / ("tmp_" + tmp_prefix)
	try:
		os.makedirs(tmp_base_folder, exist_ok=False)
		annotation_folder = tmp_base_folder / "sample_annotations"
		os.makedirs(annotation_folder, exist_ok=False)
		sample_info_with_annotations = []
		for sample_name, sample_haplotype, sample_sequence, annotation in sample_info:
			if annotation:
				Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Using annotations for sample \"{sample_name}\" haplotype \"{sample_haplotype}\" from path \"{annotation}\"", file=sys.stderr)
				sample_info_with_annotations.append((sample_name, sample_haplotype, sample_sequence, annotation))
			else:
				Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Running liftoff for sample \"{sample_name}\" haplotype \"{sample_haplotype}\"", file=sys.stderr)
				sample_annotation_file = annotation_folder / (sample_name + "_" + sample_haplotype + ".gff3.gz")
				tmp_folder = tmp_base_folder / ("tmp_" + sample_name + "_" + sample_haplotype)
				os.makedirs(tmp_folder, exist_ok=False)
				handle_new_sample_liftoff_use_tmp_folder(base_folder, tmp_folder, sample_annotation_file, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path)
				sample_info_with_annotations.append((sample_name, sample_haplotype, tmp_folder / (sample_name + "_" + sample_haplotype + ".fa"), sample_annotation_file))
		samples_with_mismatch = []
		for name, hap, fasta_file, annotation_file in sample_info_with_annotations:
			version_match = validate_gff3_is_same_isoformcheck_version(annotation_file, refannotation_hash)
			if not version_match:
				samples_with_mismatch.append((name, hap))
		if len(samples_with_mismatch) >= 1:
			if not force:
				raise Util.ParameterError("Sample annotation does not match IsoformCheck version. Rerun liftover for samples, or if you are sure about what you are doing you can force insert with --force. Samples and haplotypes with invalid versions: " + ", ".join(name + " " + hap for name, hap in samples_with_mismatch))
			else:
				print(f"{datetime.datetime.now().astimezone()}: Sample annotation does not match IsoformCheck version. Adding samples anyway due to --force. Samples and haplotypes with invalid versions: " + ", ".join(name + " " + hap for name, hap in samples_with_mismatch), file=sys.stderr)
		Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Adding sample sequences to temporary agc file", file=sys.stderr)
		temp_agc_path = add_samples_to_agc(base_folder, tmp_base_folder, sample_info_with_annotations, num_threads, agc_path)
		Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Adding sample proteins to sql database", file=sys.stderr)
		add_multiple_sample_proteins_to_database(base_folder, sample_info_with_annotations, num_threads)
		Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Copying temporary agc and annotation files to database folder", file=sys.stderr)
		copy_annotations(base_folder, sample_info_with_annotations)
		copy_agc(temp_agc_path, base_folder)
		Util.verbose_print(1, f"{datetime.datetime.now().astimezone()}: Finished adding samples successfully", file=sys.stderr)
	finally:
		shutil.rmtree(tmp_base_folder)

def copy_annotations(database_folder, sample_info_with_annotations):
	"""
	Copies temporary annotation files to the database sample annotation folder.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_info_with_annotations: List of tuples of (sample_name, sample_haplotype, sample_sequence_file, sample_annotation)
	"""
	sample_annotation_folder = database_folder / "sample_annotations"
	for sample_name, sample_haplotype, _, annotation in sample_info_with_annotations:
		cp_command = ["cp", str(annotation), str(sample_annotation_folder / (sample_name + "_" + sample_haplotype + ".gff3.gz"))]
		cp_result = subprocess.run(cp_command)
		if cp_result.returncode != 0:
			raise RuntimeError("copy did not run successfully.")

def copy_agc(temp_agc_path, database_folder):
	"""
	Copies a temporary agc file to the database agc file.

	Args:
		temp_agc_path: Path of temporary agc file
		database_folder: Folder where database files are located. Type should be pathlib.Path
	"""
	cp_command = ["cp", str(temp_agc_path), str(database_folder / "sequences.agc")]
	cp_result = subprocess.run(cp_command)
	if cp_result.returncode != 0:
		raise RuntimeError("copy did not run successfully.")

def check_if_sample_exists(database_folder, sample_name, sample_haplotype):
	"""
	Check if a sample haplotype exists.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample

	Returns:
		boolean True if sample exists, False if not
	"""
	with sqlite3.connect(str(database_folder / "sample_info.db")) as connection:
		cursor = connection.cursor()
		sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
		if sample_id is None: return False
		sample_id = sample_id[0]
		haplotype_id = cursor.execute("SELECT Id FROM Haplotype WHERE SampleId=? AND Haplotype=?", (sample_id, sample_haplotype)).fetchone()
		if haplotype_id is not None: return True
	return False

def add_samples_to_agc(database_folder, tmp_folder, sample_info, num_threads, agc_path):
	"""
	Adds multiple sample fastas to agc file, leaving a temporary agc file which must be separately moved afterwards

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		tmp_folder: Temporary folder. Final temporary agc will will be stored here. Type should be pathlib.Path
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_sequence_file, sample_annotation_file)
		num_threads: Number of threads
		agc_path: Path to agc binary

	Returns:
		Path to the completed temporary agc file in type pathlib.Path
	"""
	initial_agc_file = database_folder / "sequences.agc"
	agc_file_number = 0
	for sample_name, sample_haplotype, sample_sequence, _ in sample_info:
		agc_file_number += 1
		Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Adding sample {sample_name} haplotype {sample_haplotype} to temporary agc file", file=sys.stderr)
		next_file = tmp_folder / ("tmp_agc_" + str(agc_file_number % 2) + ".agc")
		agc_command = [agc_path, "append", "-t", str(num_threads), str(initial_agc_file), str(sample_sequence)]
		with open(str(next_file), "wb") as new_agc:
			agc_result = subprocess.run(agc_command, stdout=new_agc)
		if agc_result.returncode != 0:
			raise RuntimeError("agc did not run successfully.")
		initial_agc_file = next_file
	return next_file

def handle_new_sample_liftoff_use_tmp_folder(database_folder, tmp_folder, target_file, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path):
	"""
	Runs liftoff for a new sample and copies resulting .gff3.gz to target location.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		tmp_folder: Folder where temporary files will be stored.
		target_file: Path where the resulting .gff3.gz should be stored. Should include file ending
		sample_sequence: Sequence file of sample in fasta/fastq/gzip format
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
		liftoff_path: Path to liftoff executable
	"""

	sample_annotation_folder = database_folder / "sample_annotations"
	reference_sequence_path = database_folder / "reference.fa"
	refannotation_hash = Util.get_refannotation_hash(database_folder / "info.txt")
	temp_file_path = tmp_folder / (sample_name + "_" + sample_haplotype + ".fa")
	prepare_fasta(sample_sequence, str(temp_file_path))

	liftoff_command = [liftoff_path, "-db", str(database_folder / "reference.gff3_db"), "-o", str(tmp_folder / "tmp_annotation.gff3"), "-p", str(num_threads), "-sc", "0.95", "-copies", "-polish", "-dir", str(tmp_folder / "intermediate_files"), "-u", str(tmp_folder / "unmapped_features.txt"), str(temp_file_path), str(reference_sequence_path)]
	Util.verbose_print(1, f"Running liftoff with command:", file=sys.stderr)
	Util.verbose_print(1, f"{' '.join(liftoff_command)}", file=sys.stderr)
	liftoff_result = subprocess.run(liftoff_command, capture_output=True, text=True)
	if liftoff_result.returncode != 0:
		Util.verbose_print(0, f"Liftoff log:", file=sys.stderr)
		Util.verbose_print(0, liftoff_result.stdout, file=sys.stderr)
		Util.verbose_print(0, liftoff_result.stderr, file=sys.stderr)
		raise RuntimeError("Liftoff did not run successfully.")
	Util.verbose_print(2, f"Liftoff log:", file=sys.stderr)
	Util.verbose_print(2, liftoff_result.stdout, file=sys.stderr)
	Util.verbose_print(2, liftoff_result.stderr, file=sys.stderr)

	gzip_command = ["gzip"]
	gff3_with_version = Gff3Parser.read_gff3_as_bytes_add_isoformcheck_version_to_start(tmp_folder / "tmp_annotation.gff3_polished", refannotation_hash)
	with open(str(target_file), "wb") as compressed_gff3:
		gzip_result = subprocess.run(gzip_command, input=gff3_with_version, stdout=compressed_gff3)
		if gzip_result.returncode != 0:
			raise RuntimeError("gzip did not run successfully")

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
		FROM SampleAllele
		INNER JOIN Isoform ON Isoform.Id = SampleAllele.IsoformId
		INNER JOIN Transcript ON Transcript.Id = Isoform.TranscriptId
		INNER JOIN Gene ON Gene.Id = Transcript.GeneId
		INNER JOIN Haplotype ON Haplotype.Id = SampleAllele.HaplotypeId
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
