#!/usr/bin/env python

import gzip
import random

def make_random_prefix():
	chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
	result = ""
	for i in range(0, 10):
		result += chars[random.randint(0, len(chars)-1)]
	return result

def file_exists(filepath):
	try:
		with open(filepath) as f:
			return True
	except:
		return False

def open_maybe_gzipped(filepath):
	"""
	Open file normally or with gzip depending on extension.
	Returns a file object in text mode.
	"""
	if str(filepath).lower().endswith('.gz'):
		return gzip.open(filepath, 'rt')  # text mode
	else:
		return open(filepath, 'r')

def get_alleleset_name(alleleset):
	assert isinstance(alleleset, tuple)
	if len(alleleset) == 0: return "missing"
	return "+".join(alleleset)
