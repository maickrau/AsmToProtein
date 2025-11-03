#!/usr/bin/env python

import gzip
import random
import subprocess
import os

DBVersion = "1"

def get_git_version():
	# https://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
	# git show -s --format=%ci
	branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__))).strip().decode().strip()
	commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=os.path.dirname(os.path.abspath(__file__))).strip().decode().strip()
	time = subprocess.check_output(["git", "show", "-s", "--format=%ci"], cwd=os.path.dirname(os.path.abspath(__file__))).strip().decode().strip()
	return "branch " + branch + " commit " + commit + " " + time

Version = "development " + get_git_version()

GlobalVerbosity = 0

class ParameterError(Exception):
	pass


def get_refannotation_hash(isoformcheck_info_file):
	result = "not found"
	with open(isoformcheck_info_file) as f:
		for l in f:
			if "Filtered reference annotation md5 checksum" in l:
				result = l.strip().split(" ")[5]
	return result

def get_liftoff_version(liftoff_binary):
	liftoff_info = subprocess.run([liftoff_binary, "--version"], capture_output=True, text=True)
	liftoff_stdout = str(liftoff_info.stdout)
	return liftoff_stdout.strip()

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

def verbose_print(verbosity, s, file):
	if verbosity <= GlobalVerbosity:
		print(s, file=file)
