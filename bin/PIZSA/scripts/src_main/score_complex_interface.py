"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine calculates the score of the native structure, randonmizes the interface a 1000 times and returns the scores of the decoys structures.

INPUT - Interacting residue pairs with their corresponding weights, Pairwise potential scores
OUTPUT - Score of the native structure, Scores for 1000 decoys.
'''

"""
MODIFICATION - Used only the interface residues for scrambled decoys instead of the original set of all residues in the protein
"""

# Main module to score the native and the randomised structures
def call_score_func(chain_res, respairs, weight_dict, pot_dict):

	complex_scores = score_complex(respairs, weight_dict, pot_dict)
	score_native = complex_scores[0]
	score_respairs = complex_scores[1]
	score_decoys = get_rand_scores(chain_res, respairs, weight_dict, pot_dict)

	return score_native, score_decoys, score_respairs


# Scores the native interface
def score_complex(respairs, weight_dict, pot_dict):

	native_scores = {}
	respair_scores = {}
	#res_scores = {}

	for interface in weight_dict.keys():
		if len(weight_dict[interface]) != 0:
			score_native = 0.0
			respair_scores.update({interface : {}})
			for pair_el in weight_dict[interface]:
				resna_pair = respairs[interface][pair_el]
				weight = weight_dict[interface][pair_el]
				score = pot_dict[resna_pair] * weight
				score_native += score
				respair_scores[interface].update({pair_el : round_val(score)})

			score_native = round_val(score_native)
			native_scores.update({interface : score_native})
		else:
			pass

	# if per residue scores are needed uncomment the block below and return res_scores

	#for pair in respair_scores.keys():
	#	res_1 = pair.split('-')[0]
	#	res_2 = pair.split('-')[1]
	#	per_res_score = round_val(respair_scores[pair] / 2)

	#	if res_1 not in res_scores.keys():
	#		res_scores.update({res_1 : per_res_score})
	#	else:
	#		res_scores[res_1] += per_res_score

	#	if res_2 not in res_scores.keys():
	#		res_scores.update({res_2 : per_res_score})
	#	else:
	#		res_scores[res_2] += per_res_score

	# if per residue scores are needed add res_scores to the return statement

	return native_scores, respair_scores


# Randomises the native interface and then scores them
def get_rand_scores(chain_res, respairs, weight_dict, pot_dict):

	from random import randint

	ranres = {}
	t_ranres = {}
	score_decoys = {}

	for chain in chain_res.keys():
		ranres.update({chain : {}})
		for res in chain_res[chain].keys():
			ranres[chain].update({res : chain_res[chain][res]})

	for interface in respairs.keys():
		if len(respairs[interface]) != 0:
			score_decoys.update({interface : {}})

	for counter in range(0, 1000):
		for chain in ranres.keys():
			resnalist = ranres[chain].values()
			resnolist = ranres[chain].keys()
			for i in range(len(resnalist) - 1, -1, -1):
				j = randint(0, i)
				resnalist[i], resnalist[j] = resnalist[j], resnalist[i]
			for k in range(0, len(resnolist)):
				t_ranres.update({resnolist[k] : resnalist[k]})
		r_score_el = score_rand(t_ranres, respairs, weight_dict, pot_dict)
		for el in r_score_el:
			score_decoys[el[0]].update({counter : el[1]})

	return score_decoys


# Randomize the native interface 


# Scores the randomised interfaces
def score_rand(t_ranres, respairs, weight_dict, pot_dict):

	r_scores = []

	for interface in respairs.keys():
		if len(respairs[interface]) != 0:
			decoy_score = 0.000
			for pair in respairs[interface]:
				if pair.count('-') == 1:
					pa_res1, pa_res2 = pair.split('-')[0], pair.split('-')[1]
				else:
					split_index = pair.find('-')
					if split_index != 0:
						pa_res1, pa_res2 = pair[:split_index], pair[split_index + 1:]
					else:
						split_index = (pair.replace('-', '', 1).find('-')) + 1
						pa_res1, pa_res2 = pair[:split_index], pair[split_index + 1:]
				n_res1, n_res2 = t_ranres[pa_res1], t_ranres[pa_res2]
				n_pair = n_res1 + '-' + n_res2
				if n_res1 == 'GLY' or n_res2 == 'GLY':
					pass
				else:
					weight = weight_dict[interface][pair]
					pair_score = pot_dict[n_pair] * weight
					decoy_score += pair_score

			decoy_score = round_val(decoy_score)

			r_score_el = (interface, decoy_score)
			r_scores.append(r_score_el)

	return r_scores


# Function that rounds off a float to 3 decimal places
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num

