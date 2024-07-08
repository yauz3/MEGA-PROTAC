"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA.  If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine takes in user provided potential files (in .csv) and converts them to the appropriate format
'''

# Take in the .csv file and return the potential dict 
def convert_pot(pot_file):

	pot_dict = {}

	for line in pot_file:
		line = line.strip()
		respair = line.split()[0].replace(',', '')
		pot_val = float(line.split()[1])
		pot_val = round_val(pot_val)
		pot_dict.update({respair : pot_val})

	return pot_dict

# Rounds off the given value to 3 decimal digits
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
