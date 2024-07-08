"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine adds options to the main script for incorporating user inputs.
'''

# Provides the users flags to tweak the program
def addoptions():

	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('input_pdb', help = 'Input PDB file for the protein complex')

	parser.add_argument('-d', '--cutoff', 
						type = float,
						choices = [4.0, 6.0, 8.0],
						default = 4.0,
						help = 'Specify the distance cut-off')

	parser.add_argument('-t', '--intertype',
						choices = ["all", "mm", "ms", "ss"],
						default = "all", 
						help = 'Specify the interacting atom types')

	parser.add_argument('-o', '--outfile',
						help = 'Specify the output file')

	parser.add_argument('-cp', '--custom_pot',
						help = 'Specify the custom potential file in .csv format')

	parser.add_argument('-p1', '--protein_1',
						help = 'Specify the first protein of the interface')

	parser.add_argument('-p2', '--protein_2',
						help = 'Specify the second protein of the interface')

	parser.add_argument('-ma', '--mut_analysis',
						nargs = '*',
						default = '0',
						help = 'Specify whether Mutational Analysis is required')

	args = parser.parse_args()

	pdb_input = args.input_pdb
	cut_off = args.cutoff
	inter_type = args.intertype
	outfile = args.outfile
	protein_1 = args.protein_1
	protein_2 = args.protein_2
	custom_pot = args.custom_pot
	alscn = args.mut_analysis

	return pdb_input, cut_off, inter_type, protein_1, protein_2, outfile, custom_pot, alscn
