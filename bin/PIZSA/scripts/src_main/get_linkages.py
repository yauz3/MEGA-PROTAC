"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine parses the pairwise distance file and returns the interface residue pairs irrespective of any filters.
'''

# Get the interacting linkages for each interface
def get_linkages(dist_file):

	linkage_pairs = []
	sub1_list = set()
	sub2_list = set()

	for interface in dist_file.keys():
		for element in dist_file[interface]:
			res_1 = ':'.join(element.split()[0].split(':')[1:])
			res_1_split = res_1.split(':')
			res_1 = res_1_split[1] + ':' + res_1_split[0] + ':' + res_1_split[2]
			res_2 = ':'.join(element.split()[1].split(':')[1:])
			res_2_split = res_2.split(':')
			res_2 = res_2_split[1] + ':' + res_2_split[0] + ':' + res_2_split[2]
			#res_pair_list = sorted([res_1, res_2], key = lambda t:t[-1])
			#res_pair = res_pair_list[0] + '-' + res_pair_list[1]
			res_pair = res_1 + '-' + res_2
			sub1_list.add(res_1)
			sub2_list.add(res_2)

			if res_pair not in linkage_pairs:
				linkage_pairs.append(res_pair)

	return linkage_pairs, sub1_list, sub2_list
