"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine calculates the z-score and other metrics from the score of the native structure and the scores of the decoy structures.

INPUT - Score of the native structure & Scores of the randomized score_decoys
OUTPUT - Z-score, Avg. Background score, 
'''

# Calculates the Z-score and other reports other metrics for the multimer
def calc_zscore_mm(score_native, score_decoys):

	from math import sqrt
	import sys

	outvals = {}
	score_native_mm = 0 
	score_decoys_mm = []
	
	for interface in score_native.keys():
		score_native_mm += score_native[interface]

	for interface in score_decoys.keys():
		for element in score_decoys[interface]:
			score_decoys_mm.append(score_decoys[interface][element])

	total = sum(score_decoys_mm)
	avg_score = total / len(score_decoys_mm)

	std_dev = 0
	for score_el in score_decoys_mm:
		std_dev += (score_el - avg_score) ** 2

	std_dev = std_dev / len(score_decoys_mm)

	if std_dev == 0:
		print 'Error: Standard Deviation of 0'
		sys.exit()

	std_dev = sqrt(std_dev)

	z_score = (score_native_mm - avg_score) / std_dev

	z_score = round_val(z_score)
	avg_score = round_val(avg_score)
	std_dev = round_val(std_dev)

	outvals.update({'Multimer' : (z_score, avg_score, std_dev)})

	return score_native_mm, outvals

# Rounds off the given value to 3 decimal digits.
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
