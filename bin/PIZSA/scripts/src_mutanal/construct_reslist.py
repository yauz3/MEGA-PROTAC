"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

def construct_reslist(interface_res, mode = None):

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

 	reslist = {}

	for subunit in interface_res:
		for res in interface_res[subunit]:
			split_res = res.split(':')
			res_id = split_res[0] + ':' + interface_res[subunit][res] + ':' + split_res[1]
			if mode == 1:
				mutate_to = ['ALA']
				reslist.update({res_id	: mutate_to})
			elif mode == 2:
				mutate_to = [aa for aa in amino_acids if aa != interface_res[subunit][res]]
				reslist.update({res_id	: mutate_to})

	return reslist
