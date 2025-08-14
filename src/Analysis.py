#!/usr/bin/env python

import scipy
import DatabaseOperations

def chi_squared_test(database_file, groups, transcript):
	"""
	Performs chi-squared test between groups and allele sets.

	Args:
		database_file: Path to sql database file. Should be in pathlib format
		groups: List of groups to include
		transcript: Filter to a single transcript if given. If none, output all transcripts.

	Returns:
		Total sample size, group sample sizes and P-values per transcript. Format is (total_sample_size, [(group_name, group_sample_size)], [(transcript_name, P-value)])
	"""

	group_counts = DatabaseOperations.get_group_sample_counts(database_file, groups)
	total_sample_size = 0
	for group, count in group_counts.items():
		total_sample_size += count
	info_per_transcript = get_contingency_tables(database_file, groups, transcript)
	p_values = []
	for transcript, result in info_per_transcript:
		observed = []
		for row in result:
			for value in row[1]:
				observed.append(float(value))
		expected = []
		groups_sum = [0.0] * (len(result[0][1]))
		alleles_sum = [0.0] * len(result)
		for row in range(0, len(result)):
			assert len(result[row][1]) == len(result[0][1])
			for column in range(0, len(result[row][1])):
				groups_sum[column] += float(result[row][1][column])
				alleles_sum[row] += float(result[row][1][column])
		total_sum = 0
		for v in groups_sum:
			total_sum += v
		for row in range(0, len(result)):
			for column in range(0, len(result[row][1])):
				expected.append(alleles_sum[row] * groups_sum[column] / float(total_sum))
		test_result = scipy.stats.chisquare(f_obs=observed, f_exp=expected, ddof=len(groups_sum)+len(alleles_sum)-2)
		p_values.append((transcript, test_result.pvalue))
	return (total_sample_size, group_counts, p_values)

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
	return info_per_transcript
