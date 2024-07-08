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

# Scale cifa according to the percent of atoms that are clashing in any given interaction
# Percent of atoms clashing = atoms involved in a clash / atoms within the defined distance threshold.


import math

def sigmoid(x):
  return 1 / (1 + math.exp(-x))                            


def percent_clash(interface_residues, clashes):

	scale_dict = {}

	for interface in clashes:
		scale_dict.update({interface : {}})
		for residue_pair in clashes[interface]:
			interactors = len(interface_residues[interface][residue_pair])
			per_clash = len(clashes[interface][residue_pair]) / float(interactors)

			if per_clash < 0.6:
				cifa_scale = sigmoid(per_clash)
				scale_dict[interface].update({residue_pair : cifa_scale})
			else:
				cifa_scale = sigmoid(per_clash)
				scale_dict[interface].update({residue_pair : cifa_scale})

	return	scale_dict


def scale_cifa(weight_dict, scale_dict):

	for interface in scale_dict:
		for residue_pair in scale_dict[interface]:
			split_res = residue_pair.split('-')
			res_1 = split_res[0][4:]
			res_2 = split_res[1][4:]
			n_res_pair = res_1 + '-' + res_2
			alt_n_res_pair = res_2 + '-' + res_1

			if n_res_pair in weight_dict[interface].keys():
				old_cifa = weight_dict[interface][n_res_pair]
				weight_dict[interface].update({n_res_pair : round(float(scale_dict[interface][residue_pair]),3)})
			else:
				old_cifa = weight_dict[interface][alt_n_res_pair]
				weight_dict[interface].update({alt_n_res_pair : round(float(scale_dict[interface][residue_pair]), 3)})

	return weight_dict

