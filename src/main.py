#!/usr/bin/env python

import sys
import argparse
import pathlib
import contextlib
import math
import DatabaseOperations
import Analysis
import HandleAssembly
import Util

# https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
@contextlib.contextmanager
def open_file_or_stdout(filename=None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def run_liftoff(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-i {args.input}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-t {args.threads}", file=sys.stderr)
	print(f"--liftoff {args.liftoff}", file=sys.stderr)
	if args.output[-8:] != ".gff3.gz":
		print("-o should have file ending .gff3.gz", file=sys.stderr)
		exit(1)
	HandleAssembly.run_liftoff(pathlib.Path(args.database), args.input, args.output, int(args.threads), args.liftoff)

def add_samples(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-i {args.input}", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-t {args.threads}", file=sys.stderr)
	print(f"--liftoff {args.liftoff}", file=sys.stderr)
	print(f"--agc {args.agc}", file=sys.stderr)
	print(f"--force {args.force}", file=sys.stderr)
	HandleAssembly.handle_multiple_new_samples_liftoff_and_transcripts_from_table(pathlib.Path(args.database), args.input, int(args.threads), args.liftoff, args.agc, args.force)

def add_sample(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-i {args.input}", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-t {args.threads}", file=sys.stderr)
	print(f"--name {args.name}", file=sys.stderr)
	print(f"--haplotype {args.haplotype}", file=sys.stderr)
	print(f"--liftoff {args.liftoff}", file=sys.stderr)
	print(f"--agc {args.agc}", file=sys.stderr)
	if args.haplotype not in ['1', '2', 'mat', 'pat']:
		print("--haplotype must be one of '1', '2', 'mat', 'pat'", file=sys.stderr)
		exit(1)
	HandleAssembly.handle_one_new_sample_liftoff_and_transcripts(pathlib.Path(args.database), args.input, args.name, args.haplotype, int(args.threads), args.liftoff, args.agc)

def list_samples(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	with open_file_or_stdout(args.output) as f:
		print("Sample\tHaplotypes", file=f)
		for sample, haplotypes in DatabaseOperations.get_samples(pathlib.Path(args.database) / "sample_info.db"):
			print(f"{sample}\t{",".join(haplotypes)}", file=f)

def list_groups(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	with open_file_or_stdout(args.output) as f:
		print("Sample\tGroups", file=f)
		for sample, groups in DatabaseOperations.get_sample_groups(pathlib.Path(args.database) / "sample_info.db"):
			print(f"{sample}\t{",".join(groups)}", file=f)

def stats(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	DatabaseOperations.basic_stats(pathlib.Path(args.database) / "sample_info.db")

def add_group(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"--sample {args.sample}", file=sys.stderr)
	print(f"--group {args.group}", file=sys.stderr)
	DatabaseOperations.add_sample_to_group(pathlib.Path(args.database) / "sample_info.db", args.sample, args.group)

def remove_group(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"--sample {args.sample}", file=sys.stderr)
	print(f"--group {args.group}", file=sys.stderr)
	DatabaseOperations.remove_sample_from_group(pathlib.Path(args.database) / "sample_info.db", args.sample, args.group)

def update_names(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-r {args.database}", file=sys.stderr)
	DatabaseOperations.rename_isoforms_by_coverage(pathlib.Path(args.database) / "sample_info.db")

def contingency_table(args):
	groups = []
	if args.group:
		groups = list(args.group)
		groups.sort()
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"--transcript {args.transcript}", file=sys.stderr)
	print(f"--group {", ".join(groups)}", file=sys.stderr)
	print(f"--table {args.table}", file=sys.stderr)
	print(f"--include-gene-info {args.include_gene_info}", file=sys.stderr)
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		print("Error: Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
		exit(1)
	if len(groups) < 2 and args.table is None:
		print("Either sample table (--sample) or at least two groups (--group) are required", file=sys.stderr)
		print("Groups should all be listed at once (eg. \"--group group1 group2 group3\", not \"--group group1 --group group2 --group group3\")", file=sys.stderr)
		exit(1)
	if args.table is not None and len(groups) >= 2:
		print("Cannot use both --sample and --group parameters", file=sys.stderr)
		print("Select only either sample table (--sample) or two or more groups (--group)", file=sys.stderr)
		exit(1)
	if len(groups) >= 2:
		headergroups, info_per_transcript = Analysis.get_contingency_tables(pathlib.Path(args.database) / "sample_info.db", groups, args.transcript)
	else:
		headergroups, info_per_transcript = Analysis.get_contingency_tables_by_table(pathlib.Path(args.database) / "sample_info.db", args.table, args.transcript)
	if args.include_gene_info:
		transcript_gene_info = DatabaseOperations.get_transcript_gene_chromosome_info(pathlib.Path(args.database) / "sample_info.db")
	with open_file_or_stdout(args.output) as f:
		if args.include_gene_info:
			print(f"Chromosome\tGene\tTranscript\tAlleleset\t{"\t".join(headergroups)}", file=f)
		else:
			print(f"Transcript\tAlleleset\t{"\t".join(headergroups)}", file=f)
		for transcript, result in info_per_transcript:
			for line in result:
				if args.include_gene_info:
					print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{line[0]}\t{"\t".join(str(c) for c in line[1])}", file=f)
				else:
					print(f"{transcript}\t{line[0]}\t{"\t".join(str(c) for c in line[1])}", file=f)

def export_isoforms(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"--transcript {args.transcript}", file=sys.stderr)
	print(f"--include-gene-info {args.include_gene_info}", file=sys.stderr)
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		print("Error: Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
		exit(1)
	if args.transcript:
		result = DatabaseOperations.get_isoforms_of_transcript(pathlib.Path(args.database), args.transcript)
		result = [(args.transcript, result)]
	else:
		result = DatabaseOperations.get_isoforms_of_all_transcripts(pathlib.Path(args.database))
	if args.include_gene_info:
		transcript_gene_info = DatabaseOperations.get_transcript_gene_chromosome_info(pathlib.Path(args.database) / "sample_info.db")
	with open_file_or_stdout(args.output) as f:
		if args.include_gene_info:
			print(f"Chromosome\tGene\tTranscript\tIsoform\tCount\tSequence", file=f)
		else:
			print(f"Transcript\tIsoform\tCount\tSequence", file=f)
		for transcript, lines in result:
			for name, sequence, count in lines:
				if args.include_gene_info:
					print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{name}\t{count}\t{sequence}", file=f)
				else:
					print(f"{transcript}\t{name}\t{count}\t{sequence}", file=f)

def export_allelesets(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"--transcript {args.transcript}", file=sys.stderr)
	print(f"--include-gene-info {args.include_gene_info}", file=sys.stderr)
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		print("Error: Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
		exit(1)
	if args.transcript:
		result = DatabaseOperations.get_allelesets_of_one_transcript(pathlib.Path(args.database) / "sample_info.db", args.transcript)
	else:
		result = DatabaseOperations.get_all_allelesets(pathlib.Path(args.database) / "sample_info.db")
	if args.include_gene_info:
		transcript_gene_info = DatabaseOperations.get_transcript_gene_chromosome_info(pathlib.Path(args.database) / "sample_info.db")
	with open_file_or_stdout(args.output) as f:
		if args.include_gene_info:
			print(f"Chromosome\tGene\tTranscript\tSample\tAlleleset", file=f)
		else:
			print(f"Transcript\tSample\tAlleleset", file=f)
		for transcript, lines in result:
			for name, sequence in lines:
				if args.include_gene_info:
					print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{name}\t{sequence}", file=f)
				else:
					print(f"{transcript}\t{name}\t{sequence}", file=f)

def chi_square(args):
	groups = []
	if args.group:
		groups = list(args.group)
		groups.sort()
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"--transcript {args.transcript}", file=sys.stderr)
	print(f"--group {", ".join(groups)}", file=sys.stderr)
	print(f"--table {args.table}", file=sys.stderr)
	print(f"--include-gene-info {args.include_gene_info}", file=sys.stderr)
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		print("Error: Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
		exit(1)
	if len(groups) < 2 and args.table is None:
		print("Either sample table (--sample) or at least two groups (--group) are required", file=sys.stderr)
		print("Groups should all be listed at once (eg. \"--group group1 group2 group3\", not \"--group group1 --group group2 --group group3\")", file=sys.stderr)
		exit(1)
	if args.table is not None and len(groups) >= 2:
		print("Cannot use both --sample and --group parameters", file=sys.stderr)
		print("Select only either sample table (--sample) or two or more groups (--group)", file=sys.stderr)
		exit(1)
	if len(groups) >= 2:
		total_sample_size, group_counts, p_values = Analysis.chi_squared_test(pathlib.Path(args.database) / "sample_info.db", groups, args.transcript)
	else:
		total_sample_size, group_counts, p_values = Analysis.chi_squared_test_by_table(pathlib.Path(args.database) / "sample_info.db", args.table, args.transcript)
	print(f"Total sample size: {total_sample_size}", file=sys.stderr)
	print("Group sample sizes:", file=sys.stderr)
	for group, count in group_counts.items():
		print(f"{group}\t{count}", file=sys.stderr)
		total_sample_size += count
	print(f"Number of transcripts with samples: {len(p_values)}", file=sys.stderr)
	if args.include_gene_info:
		transcript_gene_info = DatabaseOperations.get_transcript_gene_chromosome_info(pathlib.Path(args.database) / "sample_info.db")
	with open_file_or_stdout(args.output) as f:
		if args.include_gene_info:
			print(f"Chromosome\tGene\tTranscript\tP-value", file=f)
		else:
			print(f"Transcript\tP-value", file=f)
		for transcript, p_value in p_values:
			if args.include_gene_info:
				print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{"1" if math.isnan(p_value) else str(p_value)}", file=f)
			else:
				print(f"{transcript}\t{"1" if math.isnan(p_value) else str(p_value)}", file=f)

def validate(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	errors = DatabaseOperations.check_if_haplotypes_are_fine(pathlib.Path(args.database) / "sample_info.db")
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		errors.append("Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
	if len(errors) == 0:
		print("No validation errors found, everything appears good.", file=sys.stderr)
	else:
		print(f"{len(errors)} validation errors found:", file=sys.stderr)
		for error in errors:
			print(error, file=sys.stderr)

def create_db(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-r {args.reference_genome}", file=sys.stderr)
	print(f"-a {args.annotation}", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"-t {args.threads}", file=sys.stderr)
	print(f"--liftoff {args.liftoff}", file=sys.stderr)
	print(f"--agc {args.agc}", file=sys.stderr)
	DatabaseOperations.initialize_database(pathlib.Path(args.database), args.reference_genome, args.annotation, args.liftoff, args.agc, int(args.threads))

def compare_novel(args):
	print("Input parameters:", file=sys.stderr)
	print(f"-db {args.database}", file=sys.stderr)
	print(f"--table {args.table}", file=sys.stderr)
	print(f"-t {args.threads}", file=sys.stderr)
	print(f"-o {args.output}", file=sys.stderr)
	print(f"--liftoff {args.liftoff}", file=sys.stderr)
	print(f"--force {args.force}", file=sys.stderr)
	if not DatabaseOperations.check_isoforms_have_names(args.database):
		print("Error: Some isoforms do not have names. Please run the rename command.", file=sys.stderr)
		exit(1)
	transcript_gene_info = DatabaseOperations.get_transcript_gene_chromosome_info(pathlib.Path(args.database) / "sample_info.db")
	novel_isoforms, novel_allelesets, sample_allelesets = HandleAssembly.compare_samples_to_database(pathlib.Path(args.database), args.table, int(args.threads), args.liftoff, args.force)
	with open(args.output + "_novel_isoforms.tsv", "w") as f:
		print(f"Chromosome\tGene\tTranscript\tName\tSequence", file=f)
		for transcript, name, sequence in novel_isoforms:
			print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{name}\t{sequence}", file=f)
	with open(args.output + "_novel_allelesets.tsv", "w") as f:
		print(f"Chromosome\tGene\tTranscript\tSample\tAlleleset", file=f)
		for transcript, name, alleleset in novel_allelesets:
			print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{name}\t{Util.get_alleleset_name(alleleset)}", file=f)
	with open(args.output + "_allelesets.tsv", "w") as f:
		print(f"Chromosome\tGene\tTranscript\tSample\tAlleleset", file=f)
		for transcript, name, alleleset in sample_allelesets:
			print(f"{transcript_gene_info[transcript][1]}\t{transcript_gene_info[transcript][0]}\t{transcript}\t{name}\t{Util.get_alleleset_name(alleleset)}", file=f)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="IsoformCheck", description="Protein isoform analysis from de novo genome assemblies.")
	parser.add_argument('--version', action='version', version="IsoformCheck version " + Util.Version)
	subparsers = parser.add_subparsers(dest="subparser_name")

	create_db_parser = subparsers.add_parser("initialize", description="Create new database")
	create_db_parser.add_argument('-r', '--reference-genome', required=True, help='Reference genome file (required)')
	create_db_parser.add_argument('-a', '--annotation', required=True, help='Reference annotation gff3 (required)')
	create_db_parser.add_argument('-db', '--database', required=True, help='Output database folder')
	create_db_parser.add_argument('--liftoff', default="liftoff", help="Path to liftoff")
	create_db_parser.add_argument('--agc', default="agc", help="Path to agc")
	create_db_parser.add_argument('-t', '--threads', default="4", help='Number of threads')
	create_db_parser.set_defaults(func=create_db)

	run_liftoff_parser = subparsers.add_parser("liftover", description="Lift over annotations to one haplotype")
	run_liftoff_parser.add_argument('-i', '--input', required=True, help='Haplotype sequence file (required)')
	run_liftoff_parser.add_argument('-o', '--output', required=True, help='Output annotation file')
	run_liftoff_parser.add_argument('-db', '--database', required=True, help='Database folder')
	run_liftoff_parser.add_argument('--liftoff', default="liftoff", help="Path to liftoff")
	run_liftoff_parser.add_argument('-t', '--threads', default="4", help='Number of threads')
	run_liftoff_parser.set_defaults(func=run_liftoff)

	add_sample_parser = subparsers.add_parser("addsample", description="Add a new sample")
	add_sample_parser.add_argument('-i', '--input', required=True, help='Sequence file (required)')
	add_sample_parser.add_argument('--name', required=True, help='Sample name')
	add_sample_parser.add_argument('--haplotype', required=True, help='Sample haplotype')
	add_sample_parser.add_argument('-db', '--database', required=True, help='Database folder')
	add_sample_parser.add_argument('--liftoff', default="liftoff", help="Path to liftoff")
	add_sample_parser.add_argument('--agc', default="agc", help="Path to agc")
	add_sample_parser.add_argument('-t', '--threads', default="4", help='Number of threads')
	add_sample_parser.set_defaults(func=add_sample)

	add_samples_parser = subparsers.add_parser("addsamples", description="Add multiple new samples")
	add_samples_parser.add_argument('-i', '--input', required=True, help='Sample table file (required)')
	add_samples_parser.add_argument('-db', '--database', required=True, help='Database folder')
	add_samples_parser.add_argument('--liftoff', default="liftoff", help="Path to liftoff")
	add_samples_parser.add_argument('--agc', default="agc", help="Path to agc")
	add_samples_parser.add_argument('-t', '--threads', default="4", help='Number of threads')
	add_samples_parser.add_argument('--force', action="store_true", help='Force insert samples even if validation fails')
	add_samples_parser.set_defaults(func=add_samples)

	list_sample_parser = subparsers.add_parser("listsamples", description="List all samples")
	list_sample_parser.add_argument('-db', '--database', required=True, help='Database folder')
	list_sample_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	list_sample_parser.set_defaults(func=list_samples)

	add_group_parser = subparsers.add_parser("addgroup", description="Add a sample to a group")
	add_group_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	add_group_parser.add_argument('--sample', required=True, help='Name of sample (required)')
	add_group_parser.add_argument('--group', required=True, help='Name of group (required)')
	add_group_parser.set_defaults(func=add_group)

	remove_group_parser = subparsers.add_parser("removegroup", description="Remove a sample from a group")
	remove_group_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	remove_group_parser.add_argument('--sample', required=True, help='Name of sample (required)')
	remove_group_parser.add_argument('--group', required=True, help='Name of group (required)')
	remove_group_parser.set_defaults(func=remove_group)

	list_group_parser = subparsers.add_parser("listgroups", description="List all groups per samples")
	list_group_parser.add_argument('-db', '--database', required=True, help='Database folder')
	list_group_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	list_group_parser.set_defaults(func=list_groups)

	update_names_parser = subparsers.add_parser("rename", description="Rename isoforms according to coverage")
	update_names_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	update_names_parser.set_defaults(func=update_names)

	stats_parser = subparsers.add_parser("stats", description="Print basic statistics about database")
	stats_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	stats_parser.set_defaults(func=stats)

	validate_parser = subparsers.add_parser("validate", description="Check sample haplotype validity")
	validate_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	validate_parser.set_defaults(func=validate)

	contingency_table_parser = subparsers.add_parser("contingencytable", description="Create contingency tables")
	contingency_table_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	contingency_table_parser.add_argument('--transcript', help='Name of transcript. If no transcript is given, all transcripts will be used.')
	contingency_table_parser.add_argument('--group', nargs='+', help='Names of groups (at least two required, multiple possible)')
	contingency_table_parser.add_argument('--table', help='Table with samples per group to include')
	contingency_table_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	contingency_table_parser.add_argument('--include-gene-info', action="store_true", help='Include information about gene in the output table.')
	contingency_table_parser.set_defaults(func=contingency_table)

	chi_square_parser = subparsers.add_parser("chisquare", description="Calculate chi squared P-values of group vs allele set")
	chi_square_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	chi_square_parser.add_argument('--transcript', help='Name of transcript. If no transcript is given, all transcripts will be used.')
	chi_square_parser.add_argument('--group', nargs='*', help='Names of groups to include')
	chi_square_parser.add_argument('--table', help='Table with samples per group to include')
	chi_square_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	chi_square_parser.add_argument('--include-gene-info', action="store_true", help='Include information about gene in the output table.')
	chi_square_parser.set_defaults(func=chi_square)

	export_allelesets_parser = subparsers.add_parser("exportallelesets", description="Export per-sample allele set table")
	export_allelesets_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	export_allelesets_parser.add_argument('--transcript', help='Name of transcript. If no transcript is given, all transcripts will be used.')
	export_allelesets_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	export_allelesets_parser.add_argument('--include-gene-info', action="store_true", help='Include information about gene in the output table.')
	export_allelesets_parser.set_defaults(func=export_allelesets)

	export_isoforms_parser = subparsers.add_parser("exportisoforms", description="Export isoform table")
	export_isoforms_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	export_isoforms_parser.add_argument('--transcript', help='Name of transcript. If no transcript is given, all transcripts will be used.')
	export_isoforms_parser.add_argument('-o', '--output', default="-", help='Output file (- for stdout) (default -)')
	export_isoforms_parser.add_argument('--include-gene-info', action="store_true", help='Include information about gene in the output table.')
	export_isoforms_parser.set_defaults(func=export_isoforms)

	compare_novel_parser = subparsers.add_parser("comparesamples", description="Compare samples to database")
	compare_novel_parser.add_argument('-db', '--database', required=True, help='Database folder (required)')
	compare_novel_parser.add_argument('--table', required=True, help='Table with novel samples to include')
	compare_novel_parser.add_argument('-o', '--output', default="result", help='Output prefix (default \"result\")')
	compare_novel_parser.add_argument('-t', '--threads', default="4", help='Number of threads')
	compare_novel_parser.add_argument('--liftoff', default="liftoff", help="Path to liftoff")
	compare_novel_parser.add_argument('--force', action="store_true", help='Force compare samples even if validation fails')
	compare_novel_parser.set_defaults(func=compare_novel)

	args = parser.parse_args()
	if not args.subparser_name:
		args = parser.parse_args(["--help"])
		exit(0)

	args.func(args)
