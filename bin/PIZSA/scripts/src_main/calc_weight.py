"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine calculates the weight for each residue pair.

INPUT - Interface Residue pairs with corresponding atomic distance information
OUTPUT - Weights for each interaction across the interface.
'''

# Calculates the weights for the interacting residue pairs, the weight being equal to the number of interacting atoms
# divided by the total number of atoms for that amino acid.
def get_weights(interface_respairs, options, flag_dict = None):

	inter_type = options['intertype']

	aa_all_atoms = {
		'GLY' : 4,
		'PRO' : 9, 
		'ALA' : 5, 
		'VAL' : 7,
		'LEU' : 8, 
		'ILE' : 8, 
		'MET' : 8, 
		'CYS' : 6, 
		'PHE' : 11, 
		'TYR' : 12, 
		'TRP' : 14,
		'HIS' : 10, 
		'LYS' : 9, 
		'ARG' : 11, 
		'GLN' : 9, 
		'ASN' : 8, 
		'GLU' : 9, 
		'ASP' : 8, 
		'SER' : 6, 
		'THR' : 7
		}

	aa_sc_atoms = {
		'GLY' : 0,
		'PRO' : 5, 
		'ALA' : 1, 
		'VAL' : 3,
		'LEU' : 4, 
		'ILE' : 4, 
		'MET' : 4, 
		'CYS' : 2, 
		'PHE' : 7, 
		'TYR' : 8, 
		'TRP' : 10,
		'HIS' : 6, 
		'LYS' : 5, 
		'ARG' : 7, 
		'GLN' : 5, 
		'ASN' : 4, 
		'GLU' : 5, 
		'ASP' : 4, 
		'SER' : 2, 
		'THR' : 3
		}

	aa_mc_atoms = {
		'GLY' : 4,
		'PRO' : 4, 
		'ALA' : 4, 
		'VAL' : 4,
		'LEU' : 4, 
		'ILE' : 4, 
		'MET' : 4, 
		'CYS' : 4, 
		'PHE' : 4, 
		'TYR' : 4, 
		'TRP' : 4,
		'HIS' : 4, 
		'LYS' : 4, 
		'ARG' : 4, 
		'GLN' : 4, 
		'ASN' : 4, 
		'GLU' : 4, 
		'ASP' : 4, 
		'SER' : 4, 
		'THR' : 4
		}

	main_chain_atoms = ['N', 'C', 'CA', 'O', 'OXT']

	interaction_respairs = {}

	for interface in interface_respairs.keys():
		if interface not in interaction_respairs.keys():
			interaction_respairs.update({interface : {}})
		for element in interface_respairs[interface].keys():
			if element.count('-') == 1:
				res_1, res_2 = element.split('-')[0], element.split('-')[1]
			elif element.count('-') > 1:
				res_1 = ':'.join([element.split(':')[0], element.split(':')[1], element.split(':')[2][0]]) 
				res_2 = ':'.join([element.split(':')[2][-3:], element.split(':')[3], element.split(':')[4]])
			res_1_int_atoms = set()
			res_2_int_atoms = set()
			for item in interface_respairs[interface][element]:
				res_a, res_b = ':'.join(item.split()[0].split(':')[1:]), ':'.join(item.split()[1].split(':')[1:])
				atom_a, atom_b = item.split()[0].split(':')[0], item.split()[1].split(':')[0]
				if res_a == res_1 and res_b == res_2:
					res_1_int_atoms.add(atom_a), res_2_int_atoms.add(atom_b)
				elif res_a == res_2 and res_b == res_1:
					res_1_int_atoms.add(atom_b), res_2_int_atoms.add(atom_a)
				interaction_respairs[interface].update({element : [res_1, len(res_1_int_atoms), res_2, len(res_2_int_atoms)]})

	weight_dict = {}

	for interface in interaction_respairs.keys():
		if interface not in weight_dict.keys():
			weight_dict.update({interface : {}})
		for element in interaction_respairs[interface].keys():
			if element.count('-') == 1:
				res_type_a, res_type_b = element.split('-')[0].split(':')[0], element.split('-')[1].split(':')[0]
				resno_a, resno_b = element.split('-')[0].split(':')[1], element.split('-')[1].split(':')[1]
				chain_a, chain_b = element.split('-')[0][-1], element.split('-')[1][-1]
			elif element.count('-') > 1:
				res_type_a, res_type_b = element[0:3], element.split(':')[2][-3:]
				resno_a, resno_b = element.split(':')[1], element.split(':')[3]
				chain_a, chain_b = element.split(':')[2][0], element.split(':')[-1]
			ressig_a, ressig_b = resno_a + ':' + chain_a, resno_b + ':' + chain_b
			res_1_int_atoms, res_2_int_atoms = interaction_respairs[interface][element][1], interaction_respairs[interface][element][3]

			if inter_type == 'all':
				res_1_total_atoms, res_2_total_atoms = aa_all_atoms[res_type_a], aa_all_atoms[res_type_b]
			elif inter_type == 'ss':
				res_1_total_atoms, res_2_total_atoms = aa_sc_atoms[res_type_a], aa_sc_atoms[res_type_b]
			elif inter_type == 'mm':
				res_1_total_atoms, res_2_total_atoms = aa_mc_atoms[res_type_a], aa_mc_atoms[res_type_b]
			elif inter_type == 'ms':
				flag = flag_dict[ressig_a + '-' + ressig_b]
				if flag == '1':
					res_1_total_atoms, res_2_total_atoms = aa_mc_atoms[res_type_a], aa_sc_atoms[res_type_b]
				elif flag == '2':
					res_1_total_atoms, res_2_total_atoms = aa_sc_atoms[res_type_a], aa_mc_atoms[res_type_b]

			weight = min((float(res_1_int_atoms)/float(res_1_total_atoms)), float(res_2_int_atoms)/float(res_2_total_atoms))
			weight = round_val(weight)
			weight_dict[interface].update({ressig_a + '-' + ressig_b : weight})

	return weight_dict


# Rounds off the given value to 3 decimal places
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
