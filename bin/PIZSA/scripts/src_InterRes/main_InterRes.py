"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""


def main_InterRes():

	from addoptions_InterRes import addoptions_InterRes
	from altloc_filter import altloc_check
	from cell_list_main import get_dist
	from parse_dist import parse_dist
	from get_interface_respairs import get_contacts

	options = addoptions_InterRes()

	input_dir = '/'.join(options['input_pdb'].split('/')[:-1]) + '/'
	input_file = options['input_pdb']
	input_file_extn = options['input_pdb'].split('.')[-1]

	outfile = open(input_file.replace('.'+input_file_extn, '_interespairs.txt'), 'w ')

	filtered_pdb = altloc_check(input_file)

	dist_list = get_dist(filtered_pdb)

	parsed_dist = parse_dist(dist_list)

	warnings = parsed_dist[-1]

	interResPairs = get_contacts(parsed_dist[0], 'all')
	interResPairs = interResPairs[2]

	res_seen = {}

	for interface in interResPairs:
		for respair in interResPairs[interface]:
			if respair.count('-') == 1:
				res_1 = respair.split('-')[0]
				res_2 = respair.split('-')[1]
			else:
				split_index = respair.find('-')
				if split_index != 0:
					res_1 = respair[:split_index]
					res_2 = respair[split_index + 1:]
				else:
					split_index = (respair.replace('-', '', 1).find('-')) + 1
					res_1 = respair[:split_index]
					res_2 = respair[split_index + 1:]

			res_1_sp = res_1.split(':')
			chain_1 = res_1_sp[-1]
			try:
				res_seen[chain_1]
			except KeyError:
					res_seen.update({chain_1 : []})
			res_1 = ':'.join([res_1_sp[1], res_1_sp[0], res_1_sp[2]])
			
			res_2_sp = res_2.split(':')
			chain_2 = res_2_sp[-1]
			try:
				res_seen[chain_2]
			except KeyError:
				res_seen.update({chain_2 : []})
			res_2 = ':'.join([res_2_sp[1], res_2_sp[0], res_2_sp[2]])

			if res_1 not in res_seen[chain_1]:
				res_seen[chain_1].append(res_1)

			if res_2 not in res_seen[chain_2]:
				res_seen[chain_2].append(res_2)

	for chain in sorted(res_seen.keys()):
		temp_list = sorted(res_seen[chain], key = lambda x: int(x.split(':')[0]))
		for residue in temp_list:
			outfile.write(residue+'\n')

	outfile.close()

	return
