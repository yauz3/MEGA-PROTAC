"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine checks for the occurance of alterate coordinates for atoms in the PDB.

INPUT - Atomic Coordinates in the PDB format
OUTPUT - If Alternative locations are present, Takes the coordinates with the highest occurance,
		 else if information is missing, takes the first entry.
'''

# Filters the alternative locations in the PDB file
def altloc_filter(input_lines):

	filtered_pdb = []
	altlocs = {}

	for line in input_lines:
		if line[0:4] == 'ATOM':
			line = line.strip()
			resna = line[16:20].strip()
			if len(resna) > 3:
				atom_name = line[12:16]
				chain = line[21:22]
				resnum = line[22:26]
				occupancy = line[55:61]
				sig = ':'.join([atom_name, resnum, chain])
				if sig in altlocs.keys():
					if occupancy > altlocs[sig][0]:
						line = line.replace(line[16:20], ' ' + line[17:20])
						altlocs[sig][0] = occupancy
						altlocs[sig][1] = line
					else:
						pass
				else:
					line = line.replace(line[16:20], ' ' + line[17:20])
					altlocs.update({sig : [occupancy, line]})
			else:
				filtered_pdb.append(line)

	for element in altlocs.keys():
		filtered_pdb.append(altlocs[element][1])

	return filtered_pdb


# Checks if the pdb has alternative atomic coordinates
def altloc_check(input_file):

	input_pdb = open(input_file, 'r')
	input_lines = input_pdb.readlines()

	alt_status = 0

	for line in input_lines:
		if line[0:4] == 'ATOM':
			line = line.strip()
			alt = line[16]
			if alt != ' ':
				alt_status = 1
				break
			else:
				continue

	if alt_status == 1:
		filtered_pdb = altloc_filter(input_lines)
		return filtered_pdb
	else:
		return input_lines
