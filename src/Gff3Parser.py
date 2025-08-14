#!/usr/bin/env python

import gzip
import Util

def parse_attributes(attr_string):
	attrs = {}
	for attr in attr_string.strip().split(';'):
		if '=' in attr:
			key, val = attr.split('=', 1)
			attrs[key] = val
	return attrs

def parse_gff3_gene_names_of_transcripts(gff3_path):
	"""
	Reads a GFF3 file and returns the gene names and IDs of transcripts.

	Args:
		gff3_path: Path to gff3

	Returns:
		Dictionary of transcript_id -> (gene_name, gene_ID)
	"""
	transcript_gene = {}
	gene_name = {}
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
	result = {}
	for transcript in transcript_gene:
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
	Each transcript is a dict with keys:
	  - 'transcript_id'
	  - 'protein_coding' (bool)
	  - 'exons' (list of exon dicts ordered by exon_number)
	  - 'extra_copy_number'
	  
	Each exon dict contains:
	  - 'contig', 'start', 'end', 'strand', 'exon_number'

	Returns:
	  List of transcript dicts.
	"""

	transcript_biotype = {}   # transcript_id -> biotype string
	transcript_exons = {}     # transcript_id -> list of exon dicts (with exon_number or None)
	transcript_extra_copy_number = {}

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

		transcripts.append({
			'transcript_id': tid,
			'protein_coding': protein_coding,
			'extra_copy_number': extra_copy_number,
			'exons': exons
		})

	return transcripts
