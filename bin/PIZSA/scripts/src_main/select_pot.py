"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine selects the correct potential to be used based on user inputs.

INPUT - Distance Cut-off and Type of Interaction between residues
OUTPUT - Appropriate potential given the inputs
'''

# Selects the appropriate potential
def select_pot(cut_off, inter_type):

	if cut_off == 4:
		if inter_type == 'all':
#			pot = '4_all_cmpd_cifa_avg.p'
			pot = '4_all_atomprop_matrix.p'
		elif inter_type == 'mm':
			pot = '4_mm_norm_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '4_ms_norm_cifa_avg.p'
		elif inter_type == 'ss':
#			pot = '4_ss_norm_cifa_avg.p'
			pot = '4_ss_atomprop_matrix.p'
	elif cut_off == 6:
		if inter_type == 'all':
#			pot = '6_all_cmpd_cifa_avg.p'
			pot = '6_all_atomprop_matrix.p'
		elif inter_type == 'mm':
			pot = '6_mm_cmpd_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '6_ms_cmpd_c_ij_no_avg.p'
		elif inter_type == 'ss':
#			pot = '6_ss_cmpd_cifa_no_avg.p'
			pot = '6_ss_atomprop_matrix.p'
	elif cut_off == 8:
		if inter_type == 'all':
#			pot = '8_all_norm_cifa_avg.p'
			pot = '8_all_atomprop_matrix.p'
		elif inter_type == 'mm':
			pot = '8_mm_cmpd_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '8_ms_cmpd_c_ij_no_avg.p'
		elif inter_type == 'ss':
#			pot = '8_ss_cmpd_cifa_avg.p'
			pot = '8_ss_atomprop_matrix.p'

	return pot
