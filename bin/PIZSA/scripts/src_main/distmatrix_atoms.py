"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

#Written to compute the Euclidean distances between all atoms in a pdb file

import sys
import time
import math
from numpy import array

def pdb_parse(pdb_file):
	residue = []
	chain = []
	atom_type = []
	position = []
	residue_type = []
	lines = pdb_file.readlines()
	for line in lines:
		list = line.split()
		id = list[0]
		if id == 'ATOM':
				residue.append(list[5])
				chain.append(list[4])
				atom_type.append(list[2])
				position.append(list[6:9])
				residue_type.append(list[3])
	return residue, chain, atom_type, position, residue_type

infile = sys.argv[1]

start_time = time.time()

pdb_file = open(infile)	
pos = pdb_parse(pdb_file)
atom_type, residue, chain, position, residue_type = array(pos[2]), array(pos[0]), array(pos[1]), array(pos[3]), array(pos[4])

fout = open(infile+'_distmatrix', 'w')

for x in range(0,(len(position))):
	for y in range(x, len(position)):
		if x != y:
			if residue[x] != residue[y]: 
				atom_id_A = atom_type[x] + ':' + residue_type[x] + ':' + residue[x] + ':' + chain[x]
				atom_id_B = atom_type[y] + ':' + residue_type[y] + ':' + residue[y] + ':' + chain[y]
				dist = math.sqrt(((float(position[y][0]) - float(position[x][0]))**2) + ((float(position[y][1]) - float(position[x][1]))**2) + \
				((float(position[y][2]) - float(position[x][2]))**2))
				if dist < 8:
					out = atom_id_A + '\t'+ atom_id_B + '\t' + str(dist) + '\n'
					fout.writelines(out)
fout.close()

print 'Time elapsed = ', time.time() - start_time, 's'
