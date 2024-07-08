"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the PDB file to return the number of chains etc.

INPUT - Atomic coordinates in the PDB format
OUTPUT - No. of chains, Amino acid residues in each chain, Atomic coordinates
'''

# Parses the PDB file to scan for incompatible inputs

import sys

def parse_pdb(pdb_file):

	pdb_lines = []
	chain_set = set()

	pdb = open(pdb_file, 'r').readlines()

	for line in pdb:
		if line[0:4] == 'ATOM':
			chain = line[21]
			residue_type = line[17:20]
			print residue_type
			chain_set.add(chain)
			pdb_lines.append(line)

	chain_set = sorted(chain_set)
	chain_num = len(chain_set)

	return chain_num, chain_set, pdb_lines


parsed_pdb = parse_pdb(sys.argv[1])

if len(parsed_pdb[1]) < 2:
	print "No\nInput file has just one protein"
else:
	print "Yes"

