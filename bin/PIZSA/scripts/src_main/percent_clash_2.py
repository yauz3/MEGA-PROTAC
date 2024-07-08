"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine checks for potential clashes in the query structure

INPUT - Clashes and their severity
OUTPUT - Scaled cifa values for the clashing atom types
'''

# Scale cifa according to the severity of the clash for any given interacting pairs of residue
# Clash severity defined as the maximum difference between distance threshold and interacting distances between clashing atoms.

import math 

def sigmoid(x):
	return 1 / (1 + math.exp(-x))

def percent_clash(interface_residues, clashes):

	scale_dict = {}

	for interface in clashes:
		scale_dict.update({interface : {}})
		for residue_pair in clashes[interface]:
#--- Ankit was here...
			frag1 = residue_pair.split('-')[0].split(':')[1:]
			frag2 = residue_pair.split('-')[1].split(':')[1:]
			frag1 = ':'.join(frag1)
			frag2 = ':'.join(frag2)
			clash_res_pair = '-'.join([frag1, frag2])
#--- Ankit was here...
			clash_magnitude = max(clashes[interface][residue_pair])

			scale_dict[interface].update({clash_res_pair : clash_magnitude})

	#for interface in clashes:
	#	scale_dict.update({interface : {}})
	#	for residue_pair in clashes[interface]:
	#		clash_magnitude = []
	#		for atom_clash in clashes[interface][residue_pair]:
	#			if type(atom_clash) == float:
	#				clash_severity = atom_clash
	#				clash_magnitude.append(float(clash_severity))
					#print clash_magnitude
	#		max_atom_clash = max(clash_magnitude)

	#		scale_dict[interface].update({residue_pair : max_atom_clash})

	return scale_dict


def scale_cifa(weight_dict, scale_dict):

	for interface in scale_dict:
		for residue_pair in scale_dict[interface]:
			split_residue = residue_pair.split('-')
			res_1 = split_residue[0]
			res_2 = split_residue[1]
#			print res_1, res_2
			n_res_pair = res_1 + '-' + res_2
			alt_n_res_pair = res_2 + '-' + res_1

			if n_res_pair in weight_dict[interface].keys():
#--- Ankit was here...
				old_cifa = weight_dict[interface][n_res_pair]
				if float(scale_dict[interface][residue_pair]) > 1.0:
					new_cifa = float(old_cifa) / float(scale_dict[interface][residue_pair])
				else:
					new_cifa = old_cifa
#				if new_cifa > 1.0:
#					new_cifa = 1.0
				weight_dict[interface].update({n_res_pair : new_cifa})

			else:
				old_cifa = weight_dict[interface][alt_n_res_pair]
				if float(scale_dict[interface][residue_pair]) > 1.0:
					new_cifa = float(old_cifa) / float(scale_dict[interface][residue_pair])
				else:
					new_cifa = old_cifa
#				if new_cifa > 1.0:
#					new_cifa = 1.0
				weight_dict[interface].update({alt_n_res_pair : new_cifa})
#--- Ankit was here...

	return weight_dict

