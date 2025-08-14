#!/usr/bin/env python

import sys
import argparse
import pathlib
import DatabaseOperations

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Rename alleles by coverage.")
	parser.add_argument('-db', '--database', required=True, help='Database folder (required)')

	args = parser.parse_args()

	print("Input parameters:", file=sys.stderr)
	print(f"-r {args.database}", file=sys.stderr)

	DatabaseOperations.rename_alleles_by_coverage(pathlib.Path(args.database) / "sample_info.db")
