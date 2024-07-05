"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine checks for potential clashes in the query structure

INPUT - Filtered list of interatomic distances, list of atomic van der Waal's radii
OUTPUT - Number and severity of clashes
'''

def Hbond(res_1, res_2):

	donors = ['ALA_N',
			'ARG_N', 'ARG_NE', 'ARG_NH1', 'ARG_NH2',
			'ASN_N', 'ASN_ND2',
			'ASP_N',
			'CYS_N',
			'GLN_N', 'GLN_NE2',
			'GLU_N',
			'GLY_N',
			'HIS_N', 'HIS_ND1', 'HIS_NE2',
			'ILE_N',
			'LEU_N',
			'LYS_N', 'LYS_NZ',
			'MET_N',
			'PHE_N',
			'SER_N', 'SER_OG',
			'THR_N', 'THR_OG1',
			'TRP_N', 'TRP_NE1',
			'TYR_N', 'TYR_OH',
			'VAL_N']

	acceptors = ['ALA_O',
				'ARG_O',
				'ASN_O', 'ASN_OD1',
				'ASP_O', 'ASP_OD1', 'ASP_OD2',
				'CYS_O',
				'GLN_O', 'GLN_OE1',
				'GLU_O', 'GLU_OE1', 'GLU_OE2',
				'GLY_O',
				'HIS_O', 'HIS_NE2',
				'ILE_O',
				'LEU_O',
				'LYS_O',
				'MET_O',
				'PHE_O',
				'PRO_O',
				'SER_O', 'SER_OG',
				'THR_O', 'THR_OG1',
				'TRP_O',
				'TYR_O', 'TYR_OH',
				'VAL_O']

	donate = 0
	accept = 0

	for res in (res_1, res_2):

		if res in donors:
			donate += 1
		elif res in acceptors:
			accept += 1

	if donate == 1 and accept == 1:
		return True
	else:
		return False


def find_clashes(data_dir, interface_residues, cut_off):

	import cPickle as pickle

	clashes = {}

	vdW_radii_dict = pickle.load(open(data_dir + 'vdW_radii.p', 'rb'))

	for interface in interface_residues:
		clashes.update({interface : {}})
		for res_pair in interface_residues[interface]:
			for interaction in interface_residues[interface][res_pair]:
				interaction = interaction.split()
				res_1_split = interaction[0].split(':')
				res_2_split = interaction[1].split(':')
				inter_dist = float(interaction[-1])

				res_key_1 = res_1_split[1] + '_' + res_1_split[0]
				res_key_2 = res_2_split[1] + '_' + res_2_split[0]

				hbond = Hbond(res_key_1, res_key_2)

				if hbond:
					condition = inter_dist < 0.89*(vdW_radii_dict[res_key_1] + vdW_radii_dict[res_key_2]) - 0.4

				else:
					condition = inter_dist < 0.89*(vdW_radii_dict[res_key_1] + vdW_radii_dict[res_key_2])

				# This clause to be used with percent_clash_2.py

				# Clash severity defined as the difference between the least defined distance for interactions (ie. 4 A)  
				# and the distance between clashing atoms

				if condition:
					clash_severity = 4.0 - float(inter_dist)
					if res_pair in clashes[interface].keys():
						clashes[interface][res_pair].append(clash_severity)
					else:
						clashes[interface].update({res_pair : [clash_severity]})

				# This clause to be used with percent_clash.py 

				#if condition:
				#	if res_pair in clashes[interface].keys():
				#		clashes[interface][res_pair].append(interaction)
				#	else:
				#		clashes[interface].update({res_pair : interaction})
						
	return clashes
