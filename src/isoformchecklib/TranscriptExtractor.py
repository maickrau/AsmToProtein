#!/usr/bin/env python

import os
import sys
import datetime
import isoformchecklib.Gff3Parser as Gff3Parser
import isoformchecklib.SequenceReader as SequenceReader
import isoformchecklib.Util as Util

def reverse_complement(seq):
	complement = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
	return seq.translate(complement)[::-1]

def translate_dna(dna_seq: str) -> str:
	"""
	Translate a DNA sequence into an amino acid sequence.
	
	- Stop codons are marked as 'X'.
	- Translation continues through all codons.
	- Trailing 1 or 2 leftover bases are marked as "+[acgt]" at the end of the result.
	
	Args:
		dna_seq: DNA sequence (upper/lowercase allowed)
	
	Returns:
		Amino acid sequence as string, with "+[acgt]" if trailing bases exist
	"""

	# https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
	codon_table = {
		'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
		'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
		'TTA':'L', 'TCA':'S', 'TAA':'x', 'TGA':'x',
		'TTG':'L', 'TCG':'S', 'TAG':'x', 'TGG':'W',

		'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
		'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
		'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
		'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',

		'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
		'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
		'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
		'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',

		'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
		'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
		'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
		'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G',
	}

	dna_seq = dna_seq.upper()
	n = len(dna_seq)
	full_codons = n // 3 * 3
	leftovers = n % 3

	protein_seq = []

	for i in range(0, full_codons, 3):
		codon = dna_seq[i:i+3]
		aa = codon_table.get(codon, '?')
		protein_seq.append(aa)

	protein = ''.join(protein_seq)

	if leftovers:
		protein += f"+{dna_seq[-leftovers:].lower()}"

	return protein

def get_transcript_dna_sequence(transcript, contig_sequence):
	"""
	Given a transcript dict and its contig sequence, returns the concatenated DNA sequence.

	Args:
	  transcript: dict with 'exons' key, list of exons with 'start', 'end', 'strand', 'contig'.
				  Assumes 1-based closed intervals in GFF (start and end inclusive).
	  contig_sequence: string, full contig sequence (0-based indexing).

	Returns:
	  DNA string of the spliced transcript sequence with correct strand.
	"""

	exons = transcript['exons']
	if not exons:
		return ''

	# Sort exons by exon_number to ensure correct order
	exons = sorted(exons, key=lambda e: e['exon_number'])

	# Extract sequences from contig (GFF is 1-based inclusive; Python is 0-based exclusive end)
	seq_fragments = []
	for exon in exons:
		start_0 = exon['start'] - 1  # convert to zero-based
		end_0 = exon['end']          # Python slice excludes end, so use end as-is
		exon_seq = contig_sequence[start_0:end_0]
		if exon['strand'] == "-":
			exon_seq = reverse_complement(exon_seq)
		seq_fragments.append(exon_seq)

	transcript_seq = ''.join(seq_fragments)

	return transcript_seq

def process_sample_transcripts_and_contigs(input_file, lifted_gff):
	"""
	Process a single sample to extract transcript amino acid sequences after liftoff.

	Args:
		input_file: Path to the sample fasta file.
		lifted_gff: Path to the sample lifted gff3 file.

	Returns:
		Tuple of (sample_transcripts, gene_locations, contig_lengths)
		Sample_transcripts is List[Tuple[transcript_id, transcript_aminoacid_sequence, extra_copy_number, transcript_location]]
		Gene_locations is dict of (id -> (contig, strand, start, end))
		Contig_lengths is dict of (contig -> length)
		transcript_location is tuple of (sample_contig, strand, start, end)
	"""
	import os
	import sys

	if not os.path.exists(input_file):
		raise RuntimeError(f"Sample sequence file not found")
	if not os.path.exists(lifted_gff):
		raise RuntimeError(f"Lifted annotation GFF not found")

	Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Reading annotation of sample.", file=sys.stderr)
	(transcripts, gene_locations) = Gff3Parser.parse_gff3_transcripts_with_exons(lifted_gff)
	transcript_locations = {}
	transcripts_per_contig = {}
	for transcript in transcripts:
		contig = transcript["contig"]
		transcript_locations[transcript["transcript_id"]] = (contig, transcript["strand"], transcript["start"], transcript["end"])
		if contig not in transcripts_per_contig: transcripts_per_contig[contig] = []
		transcripts_per_contig[contig].append(transcript)

	result_contig_lengths = {}
	transcripts_data = []
	Util.verbose_print(2, f"{datetime.datetime.now().astimezone()}: Reading sequences of sample.", file=sys.stderr)
	for name, contig_seq in SequenceReader.stream_sequences(input_file):
		result_contig_lengths[name] = len(contig_seq)
		if name not in transcripts_per_contig: continue
		for tx in transcripts_per_contig[name]:
			if not tx['exons']:
				continue
			seq = get_transcript_dna_sequence(tx, contig_seq)
			seq = translate_dna(seq)
			transcript_name = tx['transcript_id']
			extra_copy_number = tx.get('extra_copy_number', 0)
			if extra_copy_number != 0:
				transcript_name = "_".join(transcript_name.split("_")[:-1])
			transcripts_data.append((transcript_name, seq, extra_copy_number, transcript_locations[tx['transcript_id']]))

	return (transcripts_data, gene_locations, result_contig_lengths)
