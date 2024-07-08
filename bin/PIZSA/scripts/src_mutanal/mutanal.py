"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""


#--- Get score
def get_score(resna_pair, weight, pot_dict, iatoms):

	#--- If no atoms encountered return zero score
	if weight == 'NA':
		score = 0
		return score
	#--- Lower limit to weight
	elif weight < 0.01:
		weight = 0.01

	raw_sij = pot_dict[iatoms][resna_pair]				# interacting residue pair

	if raw_sij > 0:
		score = raw_sij * weight				# multiply weight when raw score is positive
	else:
		score = raw_sij / weight				# divide weight when raw score is negative
#	print raw_sij, weight, score

	return score



def get_interactors(mutares_id, weight_dict):

	interactors = {}

	for interface in weight_dict:
		for respair in weight_dict[interface]:
			res1 = respair.split('_')[0]
			res2 = respair.split('_')[1]

			if (mutares_id == res1) or (mutares_id == res2):
				interactors[respair] = weight_dict[interface][respair]

	return interactors



def get_score_diff(mutares_id, interactors, interface_res, mutant_el, pot_dict):

	score_diff = 0							# initialize score difference

	for respair in interactors:
		res1 = respair.split('_')[0]				# residue ID 1
		res2 = respair.split('_')[1]				# residue ID 2

		resna1 = interface_res[res1[-1]][res1]			# residue name 1
		resna2 = interface_res[res2[-1]][res2]			# residue name 2

		weight_mcmc = interactors[respair]['mcmc']		# main chain-main chain weight
		weight_mcsc = interactors[respair]['mcsc']		# main chain-side chain weight
		weight_scsc = interactors[respair]['scsc']		# side chain-side chain weight

		oripair = resna1 + '-' + resna2				# original residue pair

		ori_mcmc_score = get_score(oripair, weight_mcmc, pot_dict, 'mcmc')
		ori_mcsc_score = get_score(oripair, weight_mcsc, pot_dict, 'mcsc')
		ori_scsc_score = get_score(oripair, weight_scsc, pot_dict, 'scsc')

		ori_score = ori_mcmc_score + ori_mcsc_score + ori_scsc_score	# original score

		if res1 == mutares_id:
			mutpair = mutant_el + '-' + resna2
		else:
			mutpair = resna1 + '-' + mutant_el

		mut_mcmc_score = get_score(mutpair, weight_mcmc, pot_dict, 'mcmc')
		mut_mcsc_score = get_score(mutpair, weight_mcsc, pot_dict, 'mcsc')
		mut_scsc_score = get_score(mutpair, weight_scsc, pot_dict, 'scsc')

		mut_score = mut_mcmc_score + mut_mcsc_score + mut_scsc_score	# mutant score

		s_diff = mut_score - ori_score					# difference between mutant and original score of one interactor
		score_diff += s_diff						# add score differences of all interactors

	return score_diff



def mutanal_test(pot_dict, weight_dict, respairs_scores, interface_res, res_list, native_score, dist_cutoff):

	from collections import OrderedDict
#	from numpy import std, mean
#	import os
#	import pickle

#	parentDir = os.getcwd()
#	dataDir = parentDir + '/data/'                          		# data directory
#	background = pickle.load(open(os.path.join(dataDir, '{}ang_background.p'.format(int(dist_cutoff))), 'rb'))

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

 	score_diff_dict = OrderedDict()

 	for mutares in res_list.keys():

		mutares_split = mutares.split(':')

		mutares_id = mutares_split[0] + ':' + mutares_split[2]

		mutant_list = res_list[mutares]

		interactors = get_interactors(mutares_id, weight_dict)
#		
#
#		affctd_interfaces = []
#
#		for interface in respairs_scores:
#--- Ankit was here
# 			if mutares_id[-1] in interface:
# 			affctd_interfaces.append(interface)
#--- Ankit was here
#
 		for mutant_el in mutant_list:

			score_diff = get_score_diff(mutares_id, interactors, interface_res, mutant_el, pot_dict)

#			mutant_score  = (native_score['X-X'] + score_diff) / len(weight_dict['X-X'])
#			mutant_zscore = (mutant_score - mean(background)) / std(background)
#			print mutares, mutant_el, score_diff
#
# 			t_score_diff = {}
#
#			score_diff = 0.0
#
#			mutate_param = mutate_parameters(mutares_id, affctd_interfaces, respairs_scores, interface_res)
#			decoy_score_dict = mutate_param[0]
#			mutant_dict = mutate_param[1]
#
#			print mutares, decoy_score_dict['X-X'], mutant_dict
#
#			f_dc_score_dict = mutate_scores(decoy_score_dict, mutant_dict, mutant_el, mutares_id, pot_dict, weight_dict)
#			print mutant_el, mutares_id
#
#			for interface in f_dc_score_dict:
#				score_diff = f_dc_score_dict[interface] - native_score[interface]
#				print mutant_el, f_dc_score_dict[interface], native_score[interface], score_diff
#				t_score_diff.update({interface : score_diff})
#
#			cum_score_diff = sum(t_score_diff.values())
			score_diff_dict.update({(mutares, mutant_el) : [score_diff, {'X-X': score_diff}]})
#
#				
#		#for element in score_diff_dict:
#		#	if element[0].split(':')[1] != element[1]:
#		#		print element[0], element[1], score_diff_dict[element][0]
#
	return score_diff_dict
#
#
#
#
#def mutate_scores(decoy_score_dict, mutant_dict, mutate_to, mutares_id, pot_dict, weight_dict):
#
#	for interface in mutant_dict:
#		for element in mutant_dict[interface]:
#			pair_el = element[0]
#			partner_resna = element[1]
#--- Ankit was here...
#			m_weight_mcmc = weight_dict[interface][pair_el]['mcmc']
#			m_weight_mcsc = weight_dict[interface][pair_el]['mcsc']
#			m_weight_scsc = weight_dict[interface][pair_el]['scsc']
#
#			if pair_el.index(mutares_id) == 0:
#				m_pair = mutate_to + '-' + partner_resna
#				m_pair = mutate_to + '_' + partner_resna
#			else:
#				m_pair = partner_resna + '-' + mutate_to
#				m_pair = partner_resna + '_' + mutate_to
#
#			if 'GLY' not in m_pair:
#			m_score_mcmc = get_score(m_pair, m_weight_mcmc, pot_dict, 'mcmc')
#			m_score_mcsc = get_score(m_pair, m_weight_mcsc, pot_dict, 'mcsc')
#			m_score_scsc = get_score(m_pair, m_weight_scsc, pot_dict, 'scsc')
#
#			m_score = m_score_mcmc + m_score_mcsc + m_score_scsc
#			decoy_score_dict[interface] += m_score
#				m_score = pot_dict[m_pair] * m_weight
#				decoy_score_dict[interface] += m_score
#			else:
#				pass										# @Ankit - THIS IS IT!!!! <DANGER>
#			print pair_el, partner_resna, mutate_to
#			print m_score_mcmc, m_weight_mcmc
#			print m_score_mcsc, m_weight_mcsc
#			print m_score_scsc, m_weight_scsc
#			print ''
#	print decoy_score_dict
#--- Ankit was here...
#
#	return decoy_score_dict
#
#
#def mutate_parameters(mutares_id, affctd_interfaces, respairs_scores, interface_res):
#
#	decoy_score_dict = {}
#	mutant_dict = {}
#
#	for interface in affctd_interfaces:
#		decoy_score_dict.update({interface : 0.0 })
#--- Ankit was here...
#		print ''
#		print mutares_id
#		for pair_el in respairs_scores[interface]:
#			split_p_el = pair_el.split('-')
#			split_p_el = pair_el.split('_')
#			chain_sort = sorted([split_p_el[0][-1], split_p_el[1][-1]])
#			sp_interface = chain_sort[0] + '-' + chain_sort[1]
#			sp_interface = interface
#			if mutares_id not in pair_el:
#				previous = decoy_score_dict[sp_interface]
#				decoy_score_dict[sp_interface] += respairs_scores[sp_interface][pair_el]
#				print pair_el, decoy_score_dict[sp_interface]
#			elif pair_el.split('-')[0] == mutares_id or pair_el.split('-')[1] == mutares_id:
#			elif (pair_el.split('_')[0] == mutares_id) or (pair_el.split('_')[1] == mutares_id):
#				mutant_dict.update({sp_interface : []})						# @Ankit - THIS IS IT!!!! <DANGER>
#				if len(mutant_dict) == 0:
#					mutant_dict[sp_interface] = []
#
#				partner_res = pair_el.replace(mutares_id, '').replace('_', '')
#				partner_resna = interface_res[partner_res[-1]][partner_res]
#				mutant_dict[sp_interface].append([pair_el, partner_resna])
#				print pair_el, mutant_dict[sp_interface], partner_res, partner_resna
#			if mutares_id == '41:A':
#				print pair_el, decoy_score_dict, mutant_dict, decoy_score_dict[sp_interface] - previous
#				print '{}\t{:.3f}'.format(pair_el, decoy_score_dict[sp_interface] - previous)
#--- Ankit was here...
#
#	return decoy_score_dict, mutant_dict
