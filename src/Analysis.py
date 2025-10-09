#!/usr/bin/env python

import scipy
import DatabaseOperations

def get_groups_from_tsv(table):
	"""
	Reads a tsv file of sample groups and returns assignments of samples to groups

	Args:
		table: Path to TSV file of samples and groups

	Returns:
		Dict of (sample -> group)
	"""
	sample_in_group = {}
	row_number = 0
	with open(table) as f:
		for l in f:
			row_number += 1
			if len(l.strip()) == 0: continue
			parts = l.strip().split("\t")
			if len(parts) != 2:
				raise RuntimeError(f"Table file has wrong format. Format should be tab separated file with two columns \"Sample\" and \"Group\". Row {row} has {len(parts)} columns.")
			if row_number == 1:
				if parts[0].lower() != "sample" or parts[1].lower() != "group":
					raise RuntimeError(f"Table file has wrong format. Format should be tab separated file with two columns \"Sample\" and \"Group\". Header does not match.")
				continue
			if parts[0] in sample_in_group:
				raise RuntimeError(f"Samples should only belong to a single group. Sample {parts[0]} belongs to groups \"{sample_in_group[parts[0]]}\" and \"{parts[1]}\"")
			sample_in_group[parts[0]] = parts[1]
	return sample_in_group

def chi_squared_test_by_table(database_file, table, transcript):
	"""
	Performs chi-squared test between groups and allele sets.

	Args:
		database_file: Path to sql database file. Should be in pathlib format
		table: TSV file of samples and groups
		transcript: Filter to a single transcript if given. If none, output all transcripts.

	Returns:
		Total sample size, group sample sizes and P-values per transcript. Format is (total_sample_size, {group_name -> group_sample_size}, [(transcript_name, P-value)])
	"""
	sample_in_group = get_groups_from_tsv(table)
	temp_group_name_mapping = DatabaseOperations.add_temporary_groups(database_file, sample_in_group)
	temp_group_names = []
	result = None
	for group in temp_group_name_mapping:
		temp_group_names.append(group)
	try:
		result = chi_squared_test(database_file, temp_group_names, transcript)
	finally:
		DatabaseOperations.remove_temporary_groups(database_file, temp_group_names)
	renamed_sample_sizes = {}
	for temp_name, sample_size in result[1].items():
		renamed_sample_sizes[temp_group_name_mapping[temp_name]] = sample_size
	return (result[0], renamed_sample_sizes, result[2])

def chi_squared_test(database_file, groups, transcript):
	"""
	Performs chi-squared test between groups and allele sets.

	Args:
		database_file: Path to sql database file. Should be in pathlib format
		groups: List of groups to include
		transcript: Filter to a single transcript if given. If none, output all transcripts.

	Returns:
		Total sample size, group sample sizes and P-values per transcript. Format is (total_sample_size, {group_name -> group_sample_size}, [(transcript_name, P-value)])
	"""
	group_counts = DatabaseOperations.get_group_sample_counts(database_file, groups)
	total_sample_size = 0
	for group, count in group_counts.items():
		total_sample_size += count
	info_per_transcript = get_contingency_tables(database_file, groups, transcript)
	p_values = []
	for transcript, result in info_per_transcript[1]:
		observed = []
		for row in result:
			for value in row[1]:
				observed.append(float(value))
		expected = []
		groups_sum = [0.0] * (len(result[0][1]))
		allelesets_sum = [0.0] * len(result)
		for row in range(0, len(result)):
			assert len(result[row][1]) == len(result[0][1])
			for column in range(0, len(result[row][1])):
				groups_sum[column] += float(result[row][1][column])
				allelesets_sum[row] += float(result[row][1][column])
		total_sum = 0
		for v in groups_sum:
			total_sum += v
		for row in range(0, len(result)):
			for column in range(0, len(result[row][1])):
				expected.append(allelesets_sum[row] * groups_sum[column] / float(total_sum))
		test_result = scipy.stats.chisquare(f_obs=observed, f_exp=expected, ddof=len(groups_sum)+len(allelesets_sum)-2)
		p_values.append((transcript, test_result.pvalue))
	return (total_sample_size, group_counts, p_values)

def get_contingency_tables_by_table(database_file, table, transcript):
	"""
	Get contingency tables of group vs allele sets with groups given in a tsv file.

	Args:
		database_file: Path to sql database file. Should be in pathlib format
		table: Path to tsv file which describes groups
		transcript: Filter to a single transcript if given. If none, output all transcripts.

	Returns:
		Transcripts, their allele sets, and counts per group. Format is [(transcript, [(alleleset, [count_in_group_n])])]
	"""
	sample_in_group = get_groups_from_tsv(table)
	temp_group_name_mapping = DatabaseOperations.add_temporary_groups(database_file, sample_in_group)
	temp_group_names = []
	result = None
	for group in temp_group_name_mapping:
		temp_group_names.append(group)
	try:
		if transcript:
			transcript_result = DatabaseOperations.get_alleleset_contingency_table_by_group(database_file, transcript, temp_group_names)
			info_per_transcript = [(transcript, transcript_result)]
		else:
			info_per_transcript = DatabaseOperations.get_all_transcripts_alleleset_contingency_table_by_group(database_file, temp_group_names)
	finally:
		DatabaseOperations.remove_temporary_groups(database_file, temp_group_names)
	renamed_groups = []
	for temp_name in temp_group_names:
		renamed_groups.append(temp_group_name_mapping[temp_name])
	return renamed_groups, info_per_transcript

def get_contingency_tables(database_file, groups, transcript):
	"""
	Get contingency tables of group vs allele sets.

	Args:
		database_file: Path to sql database file. Should be in pathlib format
		groups: List of groups to include
		transcript: Filter to a single transcript if given. If none, output all transcripts.

	Returns:
		Transcripts, their allele sets, and counts per group. Format is [(transcript, [(alleleset, [count_in_group_n])])]
	"""
	if transcript:
		transcript_result = DatabaseOperations.get_alleleset_contingency_table_by_group(database_file, transcript, groups)
		info_per_transcript = [(transcript, transcript_result)]
	else:
		info_per_transcript = DatabaseOperations.get_all_transcripts_alleleset_contingency_table_by_group(database_file, groups)
	return groups, info_per_transcript
