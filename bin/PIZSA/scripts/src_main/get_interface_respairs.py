"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the pairwise distances file and returns residue pairs and the corresponding interacting atoms for each pair.

INPUT - List containing pairwise atomic distances (parsed) & Type of atomic interactions to be considered.
OUTPUT - Interface residues (chainwise), Residue pairs (interface-wise), Residue pairs with interacting atoms information
		 (to be used in calculating the weights) & Size of the Interfaces
'''

# Takes the parsed dist and the interactions_filter and returns the interacting atoms per residue pair accordingly
def get_contacts(parsed_dist, interactions_filter):

	interface_respairs = {}
	respairs = {}
	interface_res = {}
	num_interactors = {}

	main_chain_atoms = ['C', 'CA', 'O', 'OXT', 'N']

	# For interactions involving all the atoms of both the residues
	if interactions_filter == 'all':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if res_pair not in interface_respairs[interface].keys():
					interface_respairs[interface].update({res_pair : [line]})
				else:
					interface_respairs[interface][res_pair].append(line)
				resna_1, resna_2 = res_1[0:3], res_2[0:3]
				interface_res[chain_1].update({ressig_1 : resna_1})
				interface_res[chain_2].update({ressig_2 : resna_2})
				respair = ressig_1 + '-' + ressig_2
				respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the side chains of both the residues
	elif interactions_filter == 'ss':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2
			
				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})
			
				res_pair = res_1 + '-' + res_2
				if atom_1 not in main_chain_atoms and atom_2 not in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the main chain atoms of both the residues
	elif interactions_filter == 'mm':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if atom_1 in main_chain_atoms and atom_2 in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the main chain atoms of residue 1 and the side chain atoms of residue 2
	elif interactions_filter == 'ms':

		flag_dict = {}

		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if atom_1 in main_chain_atoms and atom_2 not in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})
					flag_dict.update({respair : '1'})
				elif atom_1 not in main_chain_atoms and atom_2 in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})
					flag_dict.update({respair : '2'})

	for interface in respairs.keys():
		num_interactors.update({interface : len(respairs[interface])})

	if interactions_filter != 'ms':
		return interface_res, respairs, interface_respairs, num_interactors
	else:
		return interface_res, respairs, interface_respairs, num_interactors, flag_dict
