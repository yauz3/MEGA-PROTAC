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

# Calculates the Z-score and other reports other metrics
def calc_zscore(score_native, score_decoys, output_file, dist_cutoff):

#	from math import sqrt
	from numpy import std, mean
	from err_warn_report import err_report
	import sys
	import os
	import pickle

#--- Ankit was here...
	parentDir = os.getcwd()
	dataDir = parentDir + '/data/'                          # data directory
	background = pickle.load(open(os.path.join(dataDir, '{}ang_background.p'.format(int(dist_cutoff))), 'rb'))
#--- Ankit was here...


	outvals = {}

#--- Ankit was here...
	for interface in score_decoys.keys():
#		total = 0.000
#		false_pos = 0
#		for score_el in score_decoys[interface].values():
#			total += score_el
#
#			if score_el <= score_native[interface]:
#				false_pos += 1
#
#			if score_el > score_native[interface]:
#				try:
#					min_tn
#				except NameError:
#					min_tn = score_el
#				else:
#					if score_el <= min_tn:
#						min_tn = score_el
#
#			try:
#				min_score
#			except NameError:
#				min_score = score_el
#			else:
#				if score_el < min_score:
#					min_score = score_el
#
#		false_pos_rate = false_pos / len(score_decoys[interface].values())
#		avg_score = total / len(score_decoys[interface].values())
#
#		std_dev = 0
#		for score_el in score_decoys[interface].values():
#			std_dev += (score_el - avg_score) ** 2
#	
#		std_dev = std_dev / len(score_decoys[interface].values())
#
#		if std_dev == 0:
#			error_code = 1
#			err_report(error_code, output_file)
#			sys.exit()
#
#		std_dev = sqrt(std_dev)
#
#		z_bg = []
#		for score_el in score_decoys[interface].values():
#			z_scr = (score_el - avg_score) / std_dev
#			print z_scr
#			z_bg.append(z_scr)
#
		std_dev = std(background)
		avg_score = mean(background)

		z_score = (score_native[interface] - avg_score) / std_dev
#		print str(z_score) + '\tnative'
#
#		try:
#			min_tn
#		except NameError:
#			z_min_tn = 'undef'
#			z_primer = 'undef'
#		else:
#			z_min_tn = (min_tn - avg_score) / std_dev
#			z_prime = z_score - z_min_tn
#
#		z_min = (min_score - avg_score) / std_dev
#		z_2 = z_score - z_min
#
		z_score = round_val(z_score)
		avg_score = round_val(avg_score)
		std_dev = round_val(std_dev)

		outvals.update({interface : (z_score, avg_score, std_dev)})

	return outvals

# Rounds off the given value to 3 decimal digits.
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
