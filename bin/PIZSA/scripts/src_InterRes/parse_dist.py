"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the pairwise distances file and remove unwanted entries like duplications, hydrogens etc.

INPUT - Pairwise Atomic distances file
OUTPUT - Parsed Pairwise Atomic Distances file
'''

# Parses the atomic distances file for duplications, hydrogens, nucleotide bases etc.
def parse_dist(dist_file):

	nucleotide_bases = ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT']

	amino_acids = ['GLY','PRO','ALA','VAL','LEU','ILE','MET','CYS','PHE', \
		'TYR','TRP','HIS','LYS','ARG','GLN','ASN','GLU','ASP','SER','THR']

	warnings = []

	t_parsed_dist = []
	res_seen_list = []

	parsed_dist = {}

	chain_res = {}

	for line in dist_file:
		res_1, res_2, dist = line.split()[0], line.split()[1], line.split()[2].strip()
		chain_1, chain_2 = res_1[-1], res_2[-1]

		if chain_1 not in chain_res.keys():
			chain_res.update({chain_1 : {}})

		if chain_2 not in chain_res.keys():
			chain_res.update({chain_2 : {}})

		chain_list = [chain_1, chain_2]
		chain_list.sort()

		intersig = chain_list[0] + '-' + chain_list[1]
		if intersig not in parsed_dist.keys():
			parsed_dist.update({intersig : []})

		pair_el = res_order(line)

		if pair_el not in t_parsed_dist:
			t_parsed_dist.append(pair_el)

	for el in t_parsed_dist:
		res_1, res_2, dist = el.split()[0], el.split()[1], el.split()[2]
		sorted_ch = sorted([res_1[-1], res_2[-1]])
		inter_sig = sorted_ch[0] + '-' + sorted_ch[1]
		res_pair = [res_1, res_2]

		if res_1[0] == 'H' or res_2[0] == 'H':     # Excludes Hydrogen atoms
			pass
		elif res_1[0:3] == 'OXT' or res_2[0:3] == 'OXT':    # Excludes Terminal Oxygen atoms
			pass
		elif res_pair in res_seen_list:            # Remove duplications due to alternate side chain conformations and also multiple models
			pass
		elif res_1.split(':')[1] in nucleotide_bases or res_2.split(':')[1] in nucleotide_bases:  # Remove nucleotide bases 
			warn_code = 1

			if res_1.split(':')[1] in nucleotide_bases:
				warnings.append([warn_code, res_1])
			if res_2.split(':')[1] in nucleotide_bases:
				warnings.append([warn_code, res_2])
		elif res_1.split(':')[1] not in amino_acids or res_2.split(':')[1] not in amino_acids:  # Exclude non_standard amino acids
			warn_code = 2

			if res_1.split(':')[1] not in amino_acids:
				warnings.append([warn_code, res_1])

			if res_2.split(':')[1] not in amino_acids:
				warnings.append([warn_code, res_2])
		else:
			parsed_dist[inter_sig].append(el)
			res_seen_list.append(res_pair)
			resna_1, resna_2 = res_1.split(':')[1], res_2.split(':')[1]
			ressig_1 = res_1.split(':')[2] + ':' + res_1.split(':')[3]
			ressig_2 = res_2.split(':')[2] + ':' + res_2.split(':')[3]
			chain_1, chain_2 = ressig_1[-1], ressig_2[-1]
			if ressig_1 not in chain_res[chain_1].keys():
				chain_res[chain_1].update({ressig_1 : resna_1})
			if ressig_2 not in chain_res[chain_2].keys():
				chain_res[chain_2].update({ressig_2 : resna_2})

	return parsed_dist, chain_res, warnings

# Function that defines the residue order in a residue pair
# The residues in a pair are lexicographically sorted
def res_order(line):

	res_1, res_2, dist = line.split()[0], line.split()[1], line.split()[2].strip()
	resno_1, resno_2 = res_1.split(':')[2], res_2.split(':')[2]
	chain_1, chain_2 = res_1[-1], res_2[-1]	
	ressig_1, ressig_2 = resno_1 + ':' + chain_1, resno_2 + ':' + chain_2

	if len(resno_1) == len(resno_2):
		if max(ressig_1, ressig_2) == ressig_1:
			pair_el = res_1 + '\t' + res_2 + '\t' + dist

		else:
			pair_el = res_2 + '\t' + res_1 + '\t' + dist
			
	elif len(resno_1) < len(resno_2):
		if resno_1 in resno_2:
			if resno_2.index(resno_1) == 0:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
				
			else:
				if max(ressig_1, ressig_2) == ressig_1:
					pair_el = res_1 + '\t' + res_2 + '\t' + dist
					
				else:
					pair_el = res_2 + '\t' + res_1 + '\t' + dist
					
		else:
			if max(ressig_1, ressig_2) == ressig_1:
				pair_el = res_1 + '\t' + res_2 + '\t' + dist
				
			else:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
				
	elif len(resno_2) < len(resno_1):
		if resno_2 in resno_1:
			if resno_1.index(resno_2) == 0:
				pair_el = res_1  + '\t' + res_2 + '\t' + dist
				
			else:
				if max(ressig_1, ressig_2) == ressig_1:
					pair_el = res_1 + '\t' + res_2 + '\t' + dist
					
				else:
					pair_el = res_2 + '\t' + res_1 + '\t' + dist
					
		else:
			if max(ressig_1, ressig_2) == ressig_1:
				pair_el = res_1 + '\t' + res_2 + '\t' + dist
				
			else:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
	
	return pair_el

