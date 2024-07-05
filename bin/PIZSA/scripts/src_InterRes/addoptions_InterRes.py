"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

def addoptions_InterRes():

	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('input_pdb', help = 'Input PDB file for the protein complex')

	parser.add_argument('-d', '--cutoff',
						type = float,
						choices = [4.0, 6.0, 8.0],
						default = 4.0,
						help = 'Define the distance threshold for interaction')

	args = parser.parse_args()

	pdb_input = args.input_pdb
	cut_off = args.cutoff

	options = {}

	options['input_pdb'] = pdb_input
	options['cutoff'] = cut_off

	return options
