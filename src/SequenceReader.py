#!/usr/bin/env python

import sys
import gzip
import os
import Util

def stream_sequences(filepath):
	"""
	Reads sequences from FASTA or FASTQ files, supports gzipped files.
	Returns enumerator over pairs of (name, sequence)
	"""

	with Util.open_maybe_gzipped(filepath) as f:
		first_line = f.readline()
		f.seek(0)  # rewind to start

		if first_line.startswith('>'):
			# FASTA format
			seq_lines = []
			name = None
			for line in f:
				line = line.strip()
				if line.startswith('>'):
					if name is not None:
						yield (name, "".join(seq_lines))
					name = line[1:].split()[0]
					seq_lines = []
				else:
					seq_lines.append(line)
			if name is not None:
				yield (name, "".join(seq_lines))

		elif first_line.startswith('@'):
			# FASTQ format
			while True:
				header = f.readline().strip()
				if not header:
					break  # EOF
				name = header[1:].split()[0]
				seq = f.readline().strip()
				plus = f.readline().strip()
				qual = f.readline().strip()
				yield (name, "".join(seq_lines))

		else:
			raise ValueError(f"File {filepath} is not recognized as FASTA or FASTQ.")

def read_sequences(filepath):
	"""
	Reads sequences from FASTA or FASTQ files, supports gzipped files.
	Returns list over pairs of (name, sequence)
	"""
	result = []
	for pair in stream_sequences(filepath):
		result.append(pair)
	return result
