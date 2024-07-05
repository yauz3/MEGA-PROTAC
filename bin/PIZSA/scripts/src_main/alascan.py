"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine runs the Alanine Scanning modules of the software

INPUT - Potentials, Weights for the Residue Pairs, Residues to be mutated & Native score
OUTPUT - Mutational significance of the residues wrt the potentials
'''


# Main module to run Alanine Scanning in the mode specified by the user
def call_alascan(pot_dict, weight_dict, respairs_scores, interface_res, native_score, mode = ['1'], res_list = None, dist_cutoff = '4.0'):

	import sys
	import os
	from collections import OrderedDict

	mutanal_path = os.path.abspath(os.path.join('scripts', 'src_mutanal'))
	sys.path.append(mutanal_path)

	from mutanal import mutanal_test
	from construct_reslist import construct_reslist

	if mode == ['1']:
		res_list = construct_reslist(interface_res, mode = 1)
		mutanal_results = mutanal_test(pot_dict, weight_dict, respairs_scores, interface_res, res_list, native_score, dist_cutoff)
		#alscn_results = alascan(pot_dict, weight_dict, interface_res, native_score)
	elif mode == ['2']:
		cum_effect_list = {}

		res_list = construct_reslist(interface_res, mode = 2)
		mutanal_results = mutanal_test(pot_dict, weight_dict, respairs_scores, interface_res, res_list, native_score, dist_cutoff)

		for mutation in mutanal_results:
			if mutation[0] not in cum_effect_list:
				cum_effect_list.update({mutation[0] : mutanal_results[mutation][0]})
			else:
				cum_effect_list[mutation[0]] += mutanal_results[mutation][0]

#--- Ankit was here...
#		cum_effect_dict = OrderedDict(sorted(cum_effect_list.items(), key = lambda x: x[1], reverse = True))
		cum_effect_dict = OrderedDict(sorted(cum_effect_list.items(), key = lambda x: x[1]))
#--- Ankit was here...

		sorted_mut_eff = OrderedDict()

		for el in cum_effect_dict.keys():
			temp_list = []
			for mut_el in mutanal_results:
				if mut_el[0] == el:
					temp_list.append([mut_el, mutanal_results[mut_el][0]])
#--- Ankit was here...
#			temp_list = sorted(temp_list, key = lambda x:x[-1], reverse = True)
			temp_list = sorted(temp_list, key = lambda x:x[-1])
#--- Ankit was here...
			for element in temp_list:
				sorted_mut_eff.update({element[0] : mutanal_results[element[0]]})
			
		mutanal_results = sorted_mut_eff	
		#alscn_results = alascan_all(pot_dict, weight_dict, interface_res, native_score)
	elif mode == ['3']:
		error_code = 3
		print "Error: Please provide Residue list for Mutational Analysis"
		sys.exit()
	elif len(mode) == 2 and mode[0] == '3':
		mutanal_results = mutanal_test(pot_dict, weight_dict, respairs_scores, interface_res, res_list, native_score, dist_cutoff)

	#if len(mode) == 2 and mode[0] == '3':
	#	return mutanal_results
	#else:
	#	return alscn_results

	return mutanal_results

# Mutates all the interface residues to Alanine and computes the mutational effect on the complex
def alascan(pot_dict, weight_dict, interface_res, native_score):

	decoy_score_dict = {}
	res_mutate = {}

	print interface_res

	for subunit in interface_res.keys():
		for res in interface_res[subunit].keys():
			for interface in weight_dict.keys():

				decoy_score = 0.0

				try:
					decoy_score_dict[interface]
				except KeyError:
					decoy_score_dict.update({interface : {}})

				for res_pair in weight_dict[interface].keys():
					if res in res_pair:
						res_index = res_pair.index(res)

						if res_pair.count('-') == 1:
							res_1 = res_pair.split('-')[0]
							res_2 = res_pair.split('-')[1]
						else:
							split_index = res_pair.find('-')
							if split_index != 0:
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]
							else:
								split_index = (res_pair.replace('-', '', 1).find('-')) + 1
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]

						resna_1 = interface_res[res_1[-1]][res_1]
						resna_2 = interface_res[res_2[-1]][res_2]

						if res_index == 0:
							sub_respair = 'ALA' + '-' + resna_2
						else:
							sub_respair = resna_1 + '-' + 'ALA'

						pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
						decoy_score += pair_score

					else:
						if res_pair.count('-') == 1:
							res_1 = res_pair.split('-')[0]
							res_2 = res_pair.split('-')[1]
						else:
							split_index = res_pair.find('-')
							if split_index != 0:
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]
							else:
								split_index = (res_pair.replace('-', '', 1).find('-')) + 1
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]

						resna_1 = interface_res[res_1[-1]][res_1]
						resna_2 = interface_res[res_2[-1]][res_2]

						orig_respair = resna_1 + '-' + resna_2

						pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
						decoy_score += pair_score

				decoy_score_dict[interface].update({res : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		try:
			native_score[interface]
			for el_res in decoy_score_dict[interface].keys():
				r_decoy_score = decoy_score_dict[interface][el_res]
				if r_decoy_score > native_score[interface]:
					res_mutate.update({el_res : r_decoy_score - native_score[interface]})
		except KeyError:
			pass

	return res_mutate


# Mutates all the interface residues to all the other amino acids and computes the mutational effect on the complex 
def alascan_all(pot_dict, weight_dict, interface_res, native_score):

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

 	decoy_score_dict = {}
 	res_mutate = {}

 	for subunit in interface_res.keys():
 		for res in interface_res[subunit].keys():
 			for aa in amino_acids:
 				for interface in weight_dict.keys():

 					decoy_score = 0.0

 					try:
 						decoy_score_dict[interface]
 					except KeyError:
 						decoy_score_dict.update({interface : {}})

 					try:
 						decoy_score_dict[interface][res]
 					except KeyError:
 						decoy_score_dict[interface].update({res : {}})

 					for res_pair in weight_dict[interface].keys():
 						if res in res_pair:
 							res_index = res_pair.index(res)

 							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							if res_index == 0:
								if aa != resna_1:
									sub_respair = aa + '-' + resna_2
								else:
									continue
							else:
								if aa != resna_2:
									sub_respair = resna_1 + '-' + aa
								else:
									continue

							pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
							decoy_score += pair_score

						else:
							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							orig_respair = resna_1 + '-' + resna_2

							pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
							decoy_score += pair_score

					try: 
						decoy_score_dict[interface][res][aa]
					except KeyError:
						decoy_score_dict[interface][res].update({aa : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		try:
			native_score[interface]
			for res in decoy_score_dict[interface].keys():
				score_diff = 0.0
				for aa in decoy_score_dict[interface][res].keys():
					decoy_score = decoy_score_dict[interface][res][aa]
					score_diff += decoy_score - native_score[interface]
				res_mutate.update({res : score_diff})
		except KeyError:
			pass


	return res_mutate

# OBSOLETE

# Mutates the residues submitted by the user to all other amino acids and computes their mutational effect on the complex
def mutate_reslist(pot_dict, weight_dict, interface_res, res_list, native_score):

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

 	decoy_score_dict = {}
 	res_mutate = {}

 	for subunit in res_list.keys():
 		for res in res_list[subunit].keys():
 			for aa in amino_acids:
				for interface in weight_dict.keys():

					decoy_score = 0.0

 					try:
 						decoy_score_dict[interface]
 					except KeyError:
 						decoy_score_dict.update({interface : {}})

 					try:
 						decoy_score_dict[interface][res]
 					except KeyError:
 						decoy_score_dict[interface].update({res : {}})

 					for res_pair in weight_dict[interface].keys():
 						if res in res_pair:
 							res_index = res_pair.index(res)

 							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							if res_index == 0:
								resna_1 = res_list[res_1[-1]][res_1]
								resna_2 = interface_res[res_2[-1]][res_2]
							else:
								resna_1 = interface_res[res_1[-1]][res_1]
								resna_2 = res_list[res_2[-1]][res]

							if res_index == 0:
								if aa != resna_1:
									sub_respair = aa + '-' + resna_2
							else:
								if aa != resna_2:
									sub_respair = resna_1 + '-' + aa

							pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
							decoy_score += pair_score

						else:
							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							orig_respair = resna_1 + '-' + resna_2

							pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
							decoy_score += pair_score

					try: 
						decoy_score_dict[interface][res][aa]
					except KeyError:
						decoy_score_dict[interface][res].update({aa : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		try:
			native_score[interface]
			for res in decoy_score_dict[interface].keys():
				score_diff = 0.0
				for aa in decoy_score_dict[interface][res].keys():
					decoy_score = decoy_score_dict[interface][res][aa]
					score_diff += decoy_score - native_score[interface]
				res_mutate.update({res : score_diff})
		except KeyError:
			pass


	return res_mutate


# Function that rounds off a float to 3 decimal places
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
