#!/usr/bin/env python

import os
import sys
import Gff3Parser
import SequenceReader

def reverse_complement(seq):
	complement = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
	return seq.translate(complement)[::-1]

def translate_dna(dna_seq: str) -> str:
	"""
	Translate a DNA sequence into an amino acid sequence.
	
	- Stop codons are marked as 'X'.
	- Translation continues through all codons.
	- Trailing 1 or 2 leftover bases are marked as "+1" or "+2" at the end of the result.
	
	Args:
		dna_seq: DNA sequence (upper/lowercase allowed)
	
	Returns:
		Amino acid sequence as string, with "+1" or "+2" if trailing bases exist
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
		protein += f"+{leftovers}"

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

def process_sample_transcripts(input_file, lifted_gff):
	"""
	Process a single sample to extract transcript amino acid sequences after liftoff.

	Args:
		input_file: Path to the sample fasta file.
		lifted_gff: Path to the sample lifted gff3 file.

	Returns:
		List[Tuple[transcript_id, transcript_aminoacid_sequence, extra_copy_number]]
	"""
	import os
	import sys

	if not os.path.exists(lifted_gff):
		raise RuntimeError(f"Lifted annotation GFF not found")

	print(f"{datetime.datetime.now().astimezone()}: Reading sequences of sample.", file=sys.stderr)
	contig_seqs = {name: seq for name, seq in SequenceReader.stream_sequences(input_file)}

	print(f"{datetime.datetime.now().astimezone()}: Reading annotation of sample.", file=sys.stderr)
	transcripts = Gff3Parser.parse_gff3_transcripts_with_exons(lifted_gff)

	print(f"{datetime.datetime.now().astimezone()}: Parsing transcripts.", file=sys.stderr)
	transcripts_data = []
	for tx in transcripts:
		if not tx['exons']:
			continue
		contig_name = tx['exons'][0]['contig']
		if contig_name not in contig_seqs:
			raise RuntimeError("Contig {contig_name} not found")
		contig_seq = contig_seqs[contig_name]
		seq = get_transcript_dna_sequence(tx, contig_seq)
		seq = translate_dna(seq)
		transcripts_data.append((tx['transcript_id'], seq, tx.get('extra_copy_number', 0)))

	return transcripts_data
