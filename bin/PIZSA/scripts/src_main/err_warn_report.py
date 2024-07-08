"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine writes the error and the warning files.
'''

# Write out the error files if the program aborts execution
def err_report(error_code, output_file): 

	error_file = open(output_file.replace('.out', '.err'), 'w')

	if error_code == 1:

		error_file.write("Error Code: 1 \n\nInput file does not contain a Protein Complex")

	elif error_code == 2:

		error_file.write("Error Code 2: \n\nStandard Deviation of background set is 0")

	elif error_code == 3:

		error_file.write("Error Code 3: \n\nPlease provide the residue list for Mutational Analysis")

	elif error_code == 4:

		error_file.write("Error Code 4: \n\nPlease check the format of the Mutational Analysis input file")

	elif error_code == 5:

		error_file.write("Error Code 5: \n\nResiudes provided for Mutational Analysis not found on the interface")

	elif error_code == 6:

		error_file.write("Error Code 6: \n\nNo interacting residue pairs found")

	error_file.close()

	return

# Write out the warning file if the program encounters ununsual creatures
def warn_report(warnings, output_file):

	# This functionality moveed to write_out

	return
