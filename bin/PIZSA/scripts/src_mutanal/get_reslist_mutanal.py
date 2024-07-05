"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the residue file for mutational analysis provided by the user
'''

def get_mutanal_list(filename, linkage_pairs, warnings, output_file):

	import sys
	import os

	main_path = os.path.abspath(os.path.join('scripts', 'src_main'))
	sys.path.append(main_path)

	from err_warn_report import err_report

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

	linkage_set = linkage_pairs[-2] | linkage_pairs[-1]

	fn = open(filename, 'r')

	reslist = {}

	for line in fn.readlines():
		line = line.split('\t')

		if len(line) != 2:
			error_code = 4
			err_report(error_code, output_file)
			sys.exit()

		res_mutate = line[0].strip()
		mutate_to = line[1].strip()

		if res_mutate in linkage_set:
			if len(mutate_to) == 3:
				reslist.update({res_mutate : [mutate_to]})
			elif mutate_to == '*':
				reslist.update({res_mutate : [aa for aa in amino_acids if aa != res_mutate.split(':')[1]]})
			else:
				mutate_to = mutate_to.split(',')
				for i in mutate_to: i = i.strip()
				reslist.update({res_mutate : mutate_to})
		else:
			warn_code = 3
			warnings.append([warn_code, res_mutate])

		if len(reslist) == 0:
			error_code = 5
			err_report(error_code, output_file)
			sys.exit()


	return reslist, warnings
