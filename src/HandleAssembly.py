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
import Util

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

def add_multiple_sample_proteins_to_database(database_file, sample_info):
	"""
	Adds proteins of a sample to the database. Assumes that liftoff has already been ran

	Args:
		database_file: Base folder
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_fasta_file_path, sample_annotation_file_path)
	"""
	print(f"step 1 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	all_processed_transcripts = []
	for sample_name, sample_haplotype, sample_fasta, sample_annotation in sample_info:
		all_processed_transcripts.append([])
		sample_transcripts = TranscriptExtractor.process_sample_transcripts(sample_fasta, sample_annotation)
		(gene_locations, transcript_locations) = Gff3Parser.parse_gene_transcript_locations(sample_annotation)
		for transcript_id, sequence, extra_copy in sample_transcripts:
			name = transcript_id
			if extra_copy != 0:
				name = "_".join(name.split("_")[:-1])
			all_processed_transcripts[-1].append((name, sequence, transcript_locations[transcript_id]))
	print(f"step 2 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	all_sample_contig_lens = []
	for sample_name, sample_haplotype, sample_fasta, sample_annotation in sample_info:
		all_sample_contig_lens.append({})
		for name, sequence in SequenceReader.stream_sequences(sample_fasta):
			all_sample_contig_lens[-1][name] = len(sequence)

	transcript_id_map = {}
	allele_name_map = {}
	sample_contig_db_ids = []
	all_haplotype_ids = []
	novel_alleles = set()
	with sqlite3.connect(str(database_file)) as connection:
		cursor = connection.cursor()
		cursor.arraysize = 10000
		print(f"step 3 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for db_id, transcript_id in cursor.execute("SELECT Id, Name FROM Transcript").fetchall():
			transcript_id_map[transcript_id] = db_id
		print(f"step 4 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for allele_id, transcript_id, transcript_sequence in cursor.execute("SELECT Id, TranscriptId, Sequence FROM Allele").fetchall():
			allele_name_map[(transcript_id, transcript_sequence)] = allele_id
		print(f"step 5 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		for sampleindex in range(0, len(sample_info)):
			sample_name = sample_info[sampleindex][0]
			sample_haplotype = sample_info[sampleindex][1]
			print(f"step 6 sample {sample_name} hap {sample_haplotype} time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			sample_contig_db_ids.append({})
			sample_contig_lens = all_sample_contig_lens[sampleindex]
			processed_transcripts = all_processed_transcripts[sampleindex]
			sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
			if sample_id is None:
				cursor.execute("INSERT INTO Sample (Name) VALUES (?)", (sample_name,))
				sample_id = cursor.execute("SELECT Id FROM Sample WHERE Name=?", (sample_name,)).fetchone()
				assert sample_id is not None
			sample_id = sample_id[0]
			print(f"step 7 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			cursor.execute("INSERT INTO Haplotype (SampleId, Haplotype) VALUES (?, ?)", (sample_id, sample_haplotype))
			print(f"step 8 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			haplotype_id = cursor.execute("SELECT Id FROM Haplotype WHERE SampleId=? AND Haplotype=?", (sample_id, sample_haplotype)).fetchone()
			assert haplotype_id is not None
			haplotype_id = haplotype_id[0]
			all_haplotype_ids.append(haplotype_id)
			print(f"step 9 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for name, length in sample_contig_lens.items():
				cursor.execute("INSERT INTO SampleContig (HaplotypeId, Name, Length) VALUES (?, ?, ?)", (haplotype_id, name, length))
			print(f"step 10 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for db_id, name in cursor.execute("SELECT Id, Name FROM SampleContig WHERE HaplotypeId=?", (haplotype_id,)):
				sample_contig_db_ids[sampleindex][name] = db_id
			print(f"step 11 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			for (transcript_id, sequence, location) in processed_transcripts:
				transcript_db_id = transcript_id_map[transcript_id]
				if (transcript_db_id, sequence) not in allele_name_map: novel_alleles.add((transcript_db_id, sequence))
		print(f"step 12 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		if len(novel_alleles) >= 1:
			print(f"step 13 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
			cursor.executemany("INSERT INTO Allele (TranscriptId, Sequence) VALUES (?, ?)", novel_alleles)
		print(f"step 14 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		if len(novel_alleles) >= 1:
			for allele_id, transcript_id, transcript_sequence in cursor.execute("SELECT Allele.Id, Transcript.Id, Allele.Sequence FROM Transcript INNER JOIN Allele ON Allele.TranscriptId = Transcript.Id").fetchall():
				allele_name_map[(transcript_id, transcript_sequence)] = allele_id
		print(f"step 15 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		insert_lines = []
		for sampleindex in range(0, len(sample_info)):
			processed_transcripts = all_processed_transcripts[sampleindex]
			for transcript_id, sequence, location in processed_transcripts:
				insert_lines.append((all_haplotype_ids[sampleindex], allele_name_map[(transcript_id_map[transcript_id], sequence)], sample_contig_db_ids[sampleindex][location[0]], location[1], location[2], location[3]))
		print(f"step 16 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		cursor.executemany("INSERT INTO SampleProtein (HaplotypeId, AlleleId, SampleContigId, SampleLocationStrand, SampleLocationStart, SampleLocationEnd) VALUES (?, ?, ?, ?, ?, ?)", insert_lines)
		print(f"step 17 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
		connection.commit()
		print(f"step 18 time {datetime.datetime.now().astimezone()}", file=sys.stderr)
	print(f"step 19 time {datetime.datetime.now().astimezone()}", file=sys.stderr)

def add_sample_proteins_to_database(database_file, sample_annotation, sample_fasta, sample_name, sample_haplotype):
	"""
	Adds proteins of a sample to the database. Assumes that liftoff has already been ran and sample annotation is in subfolder sample_annotation

	Args:
		database_file: Base folder
		sample_annotation: Path of sample gff3 annotation file
		sample_fasta: Path of sample fasta file
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

def handle_one_new_sample_liftoff_and_transcripts(database_folder, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path, agc_path):
	"""
	Handles one new sample. Adds the sample sequences to the agc database, runs liftoff and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_sequence: Sequence file of sample in fasta/fastq/gzip format
		sample_name: Name of new sample
		sample_haplotype: Haplotype of new sample
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
	"""
	sample_info = [(sample_name, sample_haplotype, sample_sequence, None)]
	handle_multiple_new_samples_liftoff_and_transcripts(database_folder, sample_info, num_threads, liftoff_path, agc_path)

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

def handle_multiple_new_samples_liftoff_and_transcripts_from_table(database_folder, sample_table_file, num_threads, liftoff_path, agc_path):
	"""
	Handles multiple samples. Adds the sample sequences to the agc database, runs liftoff and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.
	If sample table file has a column for annotation, then annotations are copied from there instead of rerun

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_table_file: Tsv file with descriptions of the samples
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
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
				raise RuntimeError(f"Sample table file has wrong format. File should be a tab-separated file with three columns per row. Row {row_number} has {len(parts)} columns.")
			if row_number == 1:
				if parts[0].lower() != "sample":
					raise RuntimeError(f"Sample table file has wrong format. First column should be sample name. Header instead has \"{parts[0]}\".")
				if parts[1].lower() != "haplotype":
					raise RuntimeError(f"Sample table file has wrong format. Second column should be haplotype name. Header instead has \"{parts[1]}\".")
				if parts[2].lower() != "assembly":
					raise RuntimeError(f"Sample table file has wrong format. Third column should be assembly file. Header instead has \"{parts[2]}\".")
				if len(parts) >= 4 and parts[3].lower() != "annotation":
					raise RuntimeError(f"Sample table file has wrong format. Fourth column should be annotation file. Header instead has \"{parts[3]}\".")
				has_annotation = (len(parts) == 4)
				continue
			if parts[1] not in ["1", "2", "mat", "pat"]:
				raise RuntimeError(f"Sample table file has wrong format. Row {row_number} haplotype is \"{parts[2]}\", expected one of \"1\", \"2\", \"mat\", \"pat\"")
			annotation = None
			if has_annotation:
				if len(parts) >= 4:
					if parts[3] is not None and len(parts[3]) >= 1:
						annotation = parts[3]
			sample_info.append((parts[0], parts[1], parts[2], annotation))
			if not Util.file_exists(parts[2]):
				raise RuntimeError(f"Sample assembly file in row {row_number} cannot be read: file \"{parts[2]}\", sample \"{parts[0]}\" haplotype \"{parts[1]}\"")
			if annotation is not None and not Util.file_exists(annotation):
				raise RuntimeError(f"Sample annotation file in row {row_number} cannot be read: file \"{parts[3]}\", sample \"{parts[0]}\" haplotype \"{parts[1]}\"")
	has_duplicates = False
	sample_name_and_haplotypes = set()
	for row in sample_info:
		key = (row[0], row[1])
		if key in sample_name_and_haplotypes:
			raise RuntimeError(f"Duplicate sample and haplotype: sample \"{key[0]}\" haplotype \"{key[1]}\"")
		sample_name_and_haplotypes.add(key)
	print(f"Adding {len(sample_info)} new assemblies", file=sys.stderr)
	handle_multiple_new_samples_liftoff_and_transcripts(database_folder, sample_info, num_threads, liftoff_path, agc_path)

def handle_multiple_new_samples_liftoff_and_transcripts(database_folder, sample_info, num_threads, liftoff_path, agc_path):
	"""
	Handles multiple new samples. Adds the sample sequences to the agc database, runs liftoff if necessary and gets transcripts, and adds the transcripts to the sample annotation sql database. Creates a temp folder with all temp files which is deleted in the end.

	Args:
		database_folder: Folder where database files are located. Type should be pathlib.Path
		sample_info: List of tuples of (sample_name, sample_haplotype, sample_sequence_file_path, sample_annotation_file_path)
		liftoff_path: Path to liftoff executable
		agc_path: Path to agc executable
	"""
	for sample_name, sample_haplotype, _, _ in sample_info:
		sample_exists = check_if_sample_exists(database_folder, sample_name, sample_haplotype)
		if sample_exists:
			raise RuntimeError(f"Sample already exists: name \"{sample_name}\" haplotype \"{sample_haplotype}\"")

	tmp_prefix = Util.make_random_prefix()
	tmp_base_folder = database_folder / ("tmp_" + tmp_prefix)
	try:
		os.makedirs(tmp_base_folder, exist_ok=False)
		annotation_folder = tmp_base_folder / "sample_annotations"
		os.makedirs(annotation_folder, exist_ok=False)
		sample_info_with_annotations = []
		for sample_name, sample_haplotype, sample_sequence, annotation in sample_info:
			if annotation:
				print(f"{datetime.datetime.now().astimezone()}: Using annotations for sample \"{sample_name}\" haplotype \"{sample_haplotype}\" from path \"{annotation}\"", file=sys.stderr)
				sample_info_with_annotations.append((sample_name, sample_haplotype, sample_sequence, annotation))
			else:
				print(f"{datetime.datetime.now().astimezone()}: Running liftoff for sample \"{sample_name}\" haplotype \"{sample_haplotype}\"", file=sys.stderr)
				sample_annotation_file = annotation_folder / (sample_name + "_" + sample_haplotype + ".gff3.gz")
				tmp_folder = tmp_base_folder / ("tmp_" + sample_name + "_" + sample_haplotype)
				os.makedirs(tmp_folder, exist_ok=False)
				handle_new_sample_liftoff_use_tmp_folder(database_folder, tmp_folder, sample_annotation_file, sample_sequence, sample_name, sample_haplotype, num_threads, liftoff_path)
				sample_info_with_annotations.append((sample_name, sample_haplotype, tmp_folder / (sample_name + "_" + sample_haplotype + ".fa"), sample_annotation_file))
		print(f"{datetime.datetime.now().astimezone()}: Adding sample sequences to temporary agc file", file=sys.stderr)
		temp_agc_path = add_samples_to_agc(database_folder, tmp_base_folder, sample_info_with_annotations, num_threads, agc_path)
		print(f"{datetime.datetime.now()}.astimezone(): Adding sample proteins to sql database", file=sys.stderr)
		add_multiple_sample_proteins_to_database(database_folder / "sample_info.db", sample_info_with_annotations)
		print(f"{datetime.datetime.now().astimezone()}: Copying temporary agc and annotation files to database folder", file=sys.stderr)
		copy_annotations(database_folder, sample_info_with_annotations)
		copy_agc(temp_agc_path, database_folder)
		print(f"{datetime.datetime.now().astimezone()}: Finished adding samples successfully", file=sys.stderr)
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
		print(f"{datetime.datetime.now().astimezone()}: Adding sample {sample_name} haplotype {sample_haplotype} to temporary agc file", file=sys.stderr)
		next_file = tmp_folder / ("tmp_agc_" + str(agc_file_number % 2) + ".agc")
		agc_command = [agc_path, "append", str(initial_agc_file), str(sample_sequence)]
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
	temp_file_path = tmp_folder / (sample_name + "_" + sample_haplotype + ".fa")
	prepare_fasta(sample_sequence, str(temp_file_path))

	liftoff_command = [liftoff_path, "-db", str(database_folder / "reference.gff3_db"), "-o", str(tmp_folder / "tmp_annotation.gff3"), "-p", str(num_threads), "-sc", "0.95", "-copies", "-polish", "-dir", str(tmp_folder / "intermediate_files"), "-u", str(tmp_folder / "unmapped_features.txt"), str(temp_file_path), str(reference_sequence_path)]
	print(f"Running liftoff with command:", file=sys.stderr)
	print(f"{' '.join(liftoff_command)}", file=sys.stderr)
	liftoff_result = subprocess.run(liftoff_command)
	if liftoff_result.returncode != 0:
		raise RuntimeError("Liftoff did not run successfully.")

	gzip_command = ["gzip"]
	with open(tmp_folder / "tmp_annotation.gff3") as raw_gff3:
		with open(str(target_file), "wb") as compressed_gff3:
			gzip_result = subprocess.run(gzip_command, stdin=raw_gff3, stdout=compressed_gff3)
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
