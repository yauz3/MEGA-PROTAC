"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

''' 
Get the user provided options from the command line
'''

# Parse the user provided arguments and return a dictionary 
def get_options():

	from addoptions import addoptions

	options = {}

	arguments = addoptions()

	options['input_pdb'] = arguments[0]
	options['cutoff'] = arguments[1]
	options['intertype'] = arguments[2]
	options['protein1'] = arguments[3]
	options['protein2'] = arguments[4]
	options['outfile'] = arguments[5]
	options['custom_pot'] = arguments[6]
	options['alascan'] = arguments[7]

	return options
