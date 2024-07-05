"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the residue file for alanine scanning provided by the user
'''

def get_reslist(filename, linkage_pairs, warnings):

	interface_res = set()
	linkage_set = linkage_pairs[-2] | linkage_pairs[-1]

	fn = open(filename, 'r')

	reslist = {}

	for line in fn.readlines():
		line = line.strip()
		split_res = line.split(':')
		resno = split_res[0]
		resna = split_res[1]
		chain = split_res[2]

		try:
			reslist[chain]
		except KeyError:
			reslist.update({chain : {}})

		if line in linkage_set:
			reslist[chain].update({resno + ':' + chain : resna})
		else:
			warn_code = 3
			warnings.append([warn_code, line])

	return reslist, warnings
