"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

import os
import cPickle as pickle

cur_dir = os.getcwd()

amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', 'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

GLY_list = []

for aa in amino_acids:
	if aa != 'GLY':
		pair_1 = 'GLY' + '-' + aa
		pair_2 = aa + '-' + 'GLY'
		GLY_list.append(pair_1)
		GLY_list.append(pair_2)

GLY_list.append('GLY-GLY')

for filename in os.listdir(cur_dir):
	if '.csv' in filename:
		temp_file = open(filename, 'r')
		pot_dict = {}
		for line in temp_file.readlines():
			res_pair = line.split()[0].replace(',', '')
			score = round(float(line.split()[1]),3)
			pot_dict.update({res_pair : score})
		if 'ss' in filename:
			for element in GLY_list:
				pot_dict.update({element : 0.000})
		pot_filename = filename.replace('.csv', '.p')
		pot_filename = pot_filename.replace('_scores', '')
		pickle.dump(pot_dict, open(pot_filename, 'wb'))
		temp_file.close()
