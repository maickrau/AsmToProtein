#!/usr/bin/env python

import gzip
import datetime
import Util

def parse_attributes(attr_string):
	attrs = {}
	for attr in attr_string.strip().split(';'):
		if '=' in attr:
			key, val = attr.split('=', 1)
			attrs[key] = val
	return attrs

def read_gff3_as_bytes_add_isoformcheck_version_to_start(input_gff3_path, filtered_refannotation_hash):
	"""
	Reads a gff3 file as bytes. Add IsoformCheck DB version and hash of annotation to the end of the comments at start of file.

	Args:
		input_gff3_path: Path to input gff3
		filtered_refannotation_hash: String with hash of annotation file

	Returns:
		Bytes which represents the file plus possibly the line about IsoformCheck version
	"""
	result_lines = []
	added_information = False
	with Util.open_maybe_gzipped(input_gff3_path) as f:
		for l in f:
			if l.startswith('#'):
				result_lines.append(l.strip())
			else:
				if not added_information:
					result_lines.append("# IsoformCheck DB version " + Util.DBVersion + " annotation hash " + filtered_refannotation_hash)
					added_information = True
				result_lines.append(l.strip())
	big_string = "\n".join(result_lines)
	result = big_string.encode('utf-8')
	return result

def filter_gff3_to_things_with_CDS(input_gff3_path, output_gff3_path):
	"""
	Filters a gff3 file to genes, transcripts, and CDSses of only the transcripts which have CDS

	Args:
		input_gff3_path: Path to input gff3
		output_gff3_path: Path to output gff3
	"""
	transcript_gene = parse_gff3_gene_names_of_transcripts(input_gff3_path)
	valid_transcripts = set()
	valid_genes = set()
	for transcript, (gene_name, gene_id) in transcript_gene.items():
		valid_transcripts.add(transcript)
		valid_genes.add(gene_id)
	header_printed = False
	with Util.open_maybe_gzipped(input_gff3_path) as f:
		with open(output_gff3_path, "w") as out_f:
			for line in f:
				if line.startswith('#'):
					out_f.write(line)
					continue
				if not header_printed:
					out_f.write(f"# IsoformCheck filter {datetime.datetime.now().astimezone()} DB version {Util.DBVersion}\n")
					header_printed = True
				if not line.strip():
					continue
				fields = line.strip().split('\t')
				if len(fields) != 9:
					continue
				seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
				if feature_type not in ["transcript", "gene", "CDS"]:
					continue
				attrs = parse_attributes(attributes)
				transcript_id = None
				gene_id = None
				if feature_type == "transcript":
					transcript_id = attrs['ID']
					gene_id = attrs['Parent']
				if feature_type == "gene":
					gene_id = attrs['ID']
				if feature_type == "CDS":
					transcript_id = attrs['Parent']
				if transcript_id in valid_transcripts or gene_id in valid_genes:
					out_f.write(line)

def parse_gff3_gene_names_of_transcripts(gff3_path):
	"""
	Reads a GFF3 file and returns the gene names and IDs of transcripts which have coding sequences.

	Args:
		gff3_path: Path to gff3

	Returns:
		Dictionary of transcript_id -> (gene_name, gene_ID)
	"""
	transcript_gene = {}
	gene_name = {}
	transcripts_with_CDS = set()
	with Util.open_maybe_gzipped(gff3_path) as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			fields = line.strip().split('\t')
			if len(fields) != 9:
				continue
			seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
			attrs = parse_attributes(attributes)
			if feature_type == "transcript":
				transcript_gene[attrs['ID']] = attrs['Parent']
			if feature_type == "gene":
				gene_name[attrs['ID']] = attrs['gene_name']
			if feature_type == "CDS":
				transcripts_with_CDS.add(attrs['Parent'])

	result = {}
	for transcript in transcript_gene:
		if transcript not in transcripts_with_CDS: continue
		result[transcript] = (gene_name[transcript_gene[transcript]], transcript_gene[transcript])
	return result

def parse_gene_transcript_locations(gff3_path):
	"""
	Parses gene and transcript locations from GFF3.

	Args:
		gff3_path: Gff3 path

	Returns:
		(gene_locations, transcript_locations) where both gene_locations and transcript_locations are dicts of (id -> (contig, strand, start, end))
	"""
	gene_locations = {}
	transcript_locations = {}
	with Util.open_maybe_gzipped(gff3_path) as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			fields = line.strip().split('\t')
			if len(fields) != 9:
				continue
			seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
			attrs = parse_attributes(attributes)
			if feature_type == "transcript":
				transcript_locations[attrs['ID']] = (seqid, strand, int(start), int(end))
			if feature_type == "gene":
				gene_locations[attrs['ID']] = (seqid, strand, int(start), int(end))

	return (gene_locations, transcript_locations)

def parse_gff3_transcripts_with_exons(gff3_path):
	"""
	Parses a GFF3 file and returns a list of transcripts.

	Args:
		gff3_path: Path to gff3(.gz) file

	Returns:
		(transcript_info, gene_info)
		Transcript_info is a dict with keys:
		  - 'transcript_id'
		  - contig
		  - start
		  - end
		  - 'protein_coding' (bool)
		  - 'exons' (list of exon dicts ordered by exon_number)
		  - 'extra_copy_number'
		Each exon dict contains:
		  - 'contig', 'start', 'end', 'strand', 'exon_number'
		Gene_info is a dict of (gene -> (contig, start, end))
	"""

	transcript_biotype = {}   # transcript_id -> biotype string
	transcript_exons = {}     # transcript_id -> list of exon dicts (with exon_number or None)
	transcript_extra_copy_number = {}
	transcript_locations = {}
	gene_locations = {}

	with Util.open_maybe_gzipped(gff3_path) as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			fields = line.strip().split('\t')
			if len(fields) != 9:
				continue

			seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
			attrs = parse_attributes(attributes)

			if feature_type in ['transcript', 'mRNA', 'gene']:
				tid = attrs.get('ID') or attrs.get('transcript_id')
				extra_copy_number = attrs.get('extra_copy_number', "0")
				biotype = attrs.get('gene_biotype') or attrs.get('biotype') or attrs.get('gene_type') or ''
				if tid:
					transcript_biotype[tid] = biotype.lower()
					transcript_extra_copy_number[tid] = int(extra_copy_number)
			if feature_type == "transcript":
				tid = attrs.get('ID') or attrs.get('transcript_id')
				transcript_locations[tid] = (seqid, strand, int(start), int(end))
			if feature_type == "gene":
				gene_id = attrs.get('ID')
				gene_locations[gene_id] = (seqid, strand, int(start), int(end))

			if feature_type == 'CDS':
				parent = attrs.get('Parent')
				if not parent:
					continue
				transcript_id = parent.split(',')[0]
				exon_number = attrs.get('exon_number')
				try:
					exon_number = int(exon_number) if exon_number is not None else None
				except ValueError:
					exon_number = None

				exon_info = {
					'contig': seqid,
					'start': int(start),
					'end': int(end),
					'strand': strand,
					'exon_number': exon_number
				}

				transcript_exons.setdefault(transcript_id, []).append(exon_info)

	transcripts = []
	for tid, exons in transcript_exons.items():
		# Sort exons: by exon_number if present, else by start coordinate
		if any(exon['exon_number'] is None for exon in exons):
			exons = sorted(exons, key=lambda e: e['start'])
			# Assign exon_number if missing
			for i, exon in enumerate(exons, 1):
				exon['exon_number'] = exon['exon_number'] or i
		else:
			exons = sorted(exons, key=lambda e: e['exon_number'])

		protein_coding = (transcript_biotype.get(tid, '') == 'protein_coding')
		extra_copy_number = transcript_extra_copy_number[tid]
		contig = transcript_locations[tid][0]
		strand = transcript_locations[tid][1]
		start = transcript_locations[tid][2]
		end = transcript_locations[tid][3]

		transcripts.append({
			'transcript_id': tid,
			'contig': contig,
			'strand': strand,
			'start': start,
			'end': end,
			'protein_coding': protein_coding,
			'extra_copy_number': extra_copy_number,
			'exons': exons
		})

	return (transcripts, gene_locations)
