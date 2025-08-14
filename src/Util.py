#!/usr/bin/env python

import gzip

def open_maybe_gzipped(filepath):
	"""
	Open file normally or with gzip depending on extension.
	Returns a file object in text mode.
	"""
	if str(filepath).lower().endswith('.gz'):
		return gzip.open(filepath, 'rt')  # text mode
	else:
		return open(filepath, 'r')
