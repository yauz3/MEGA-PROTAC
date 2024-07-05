"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This script contains the final predict_binding routines.

INPUT - Atomic Coordinates in PDB format
OUTPUT - Scores file, Interacting Residues file and Alanine Scanning results file
'''

# Parse dist file in case the software is run for all the interfaces of a multimeric complex
def parse_dist_all(pdb_file, cut_off):

	from altloc_filter import altloc_check
	from cell_list_main import get_dist
	from parse_dist import parse_dist

	filtered_pdb = altloc_check(pdb_file)

	dist_list = get_dist(filtered_pdb, cut_off)

	parsed_dist = parse_dist(dist_list)

	return parsed_dist

# Main module to run the software on all the interfaces of a multimeric complex
def predict_binding_all(parent_dir, input_file, pot_dict, z_threshold, output_file, options, f_alascan = None):

	import os
	import sys 
	from get_interface_respairs import get_contacts
	from find_clashes_3 import find_clashes
	from percent_clash_3 import percent_clash, scale_cifa
	from calc_weight import get_weights
	from get_linkages import get_linkages
#	from score_complex import call_score_func
	from score_complex_clash import call_score_func
	from calc_zscore import calc_zscore
	from calc_zscore_mm import calc_zscore_mm
	from alascan import call_alascan
	from get_reslist_als import get_reslist
	from err_warn_report import err_report
	from write_out import write_out
	import get_interaction_propensity as intprop

	mutanal_path = os.path.abspath(os.path.join('scripts', 'src_mutanal'))
	sys.path.append(mutanal_path)

	from get_reslist_mutanal import get_mutanal_list
	from mutanal import mutanal_test

	pdb_name = output_file.split('/')[-1].replace('.out', '')

	parsed_dist = parse_dist_all(input_file, options['cutoff'])

	warnings = parsed_dist[-1]

	interface_residue_pairs = get_contacts(parsed_dist[0], options['intertype'])
#	print parsed_dist[0]
#	interface_prop = intprop.getPropensity(parsed_dist[0])
	interface_prop, opt_alpha_fraction = intprop.getPropensity(parsed_dist[0], options['cutoff'], options['intertype'])

# Ankit was here...
	num_interactors = interface_residue_pairs[3]

	count = 0
	for keys in num_interactors.keys():
		count += num_interactors[keys]

	mod_num_interactors = {}
	mod_num_interactors['X-X'] = count
# Ankit was here...

	if len(num_interactors) == 0:
		error_code = 6
		err_report(error_code, output_file)
		sys.exit()

# Ankit was here...
	count = 0

	premod_2 = {}

	for keys in interface_residue_pairs[2].keys():
		count += 1
#		print "=============== {} ==============\n\n\n".format(count)
#		print keys
		premod_2.update(interface_residue_pairs[2][keys])

	mod_2 = {}
	mod_2['X-X'] = premod_2
#	print mod_2
#	intprop.getPropensity(mod_2)

#--- Ankit again!
#--- Adding clash detection here...
	clashes = find_clashes(parent_dir + '/data/', mod_2, options['cutoff'])

	scale_dict = percent_clash(mod_2, clashes)
#--- Clashes detected!!!


# Ankit was here...

	if options['intertype'] == 'ms':
		weight_dict = get_weights(interface_residue_pairs[2], options, flag_dict = interface_residue_pairs[4])
	else:
		weight_dict = get_weights(interface_residue_pairs[2], options)

# Ankit was here...
		mod_weight_dict = get_weights(mod_2, options)

#--- Ankit again!
#--- Rescaling weights here...
#	mod_s_weight_dict = scale_cifa(mod_weight_dict, scale_dict)
	mod_s_weight_dict = scale_cifa(mod_weight_dict, scale_dict, interface_prop)
#--- Weights rescaled!

	premod_1 = {}
	for keys in interface_residue_pairs[1].keys():
		premod_1.update(interface_residue_pairs[1][keys])

	mod_1 = {}
	mod_1['X-X'] = premod_1


	premod_parse = {}
	for keys in parsed_dist[1].keys():
		premod_parse.update(parsed_dist[1][keys])

	mod_parse = {}
	mod_parse['X'] = premod_parse
# Ankit was here...

	#scores = call_score_func(parsed_dist[1], interface_residue_pairs[1], weight_dict, pot_dict)
# Ankit was here...
	mod_scores = call_score_func(mod_parse, mod_1, mod_s_weight_dict, pot_dict, scale_dict)
# Ankit was here...

	#native_score = scores[0]
	#decoy_scores = scores[1]
	#respair_scores = scores[2]

# Ankit was here...
	mod_native_score = mod_scores[0]
	mod_decoy_scores = mod_scores[1]
	mod_respair_scores = mod_scores[2]
	int_len = len(mod_respair_scores['X-X'])
	norm_score = {'X-X' : mod_native_score['X-X'] / float(int_len)}
# Ankit was here...

	#print native_score.keys()
	#print decoy_scores.keys()

#	out_val = calc_zscore(native_score, decoy_scores, output_file)
#	print out_val
	#out_val_mm = calc_zscore_mm(native_score, decoy_scores)

# Ankit was here...
	mod_out_val = calc_zscore(norm_score, mod_decoy_scores, output_file, options['cutoff'])
#	print mod_out_val
# Ankit was here...

	#native_score.update({'Multimer' : out_val_mm[0]})
	#out_val.update({'Multimer' : out_val_mm[1]})

	linkage_pairs = get_linkages(parsed_dist[0])

	# Replace alascan with mutanal later

	if f_alascan == None:
		pass
	elif len(f_alascan) == 2:
		reslist = get_mutanal_list(f_alascan[1], linkage_pairs, warnings, output_file)
		#mutanal = mutanal_test(pot_dict, weight_dict, respair_scores, interface_residue_pairs[0], reslist[0], native_score)
		#reslist = get_reslist(f_alascan[1], linkage_pairs, warnings)
		warnings = reslist[-1]
		alascan = call_alascan(pot_dict, weight_dict, respair_scores, interface_residue_pairs[0], native_score, mode = f_alascan, res_list = reslist[0], dist_cutoff = options['cutoff'])
	else:
#--- Ankit was here...
#		alascan = call_alascan(pot_dict, weight_dict, respair_scores, interface_residue_pairs[0], native_score, mode = f_alascan)
		alascan = call_alascan(pot_dict, mod_s_weight_dict, mod_respair_scores, interface_residue_pairs[0], mod_native_score, mode = f_alascan, dist_cutoff = options['cutoff'])
#--- Ankit was here...

	if f_alascan == None:
		#output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], respair_scores, warnings)
		mod_output = write_out(output_file, mod_out_val, options, mod_num_interactors, mod_native_score, z_threshold, mod_1, mod_respair_scores, warnings, opt_alpha_fraction)
	else:
#		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], respair_scores, warnings, alscn = alascan)
		mod_output = write_out(output_file, mod_out_val, options, mod_num_interactors, mod_native_score, z_threshold, mod_1, mod_respair_scores, warnings, opt_alpha_fraction, alscn = alascan)

	return 

# Main module to run the software for dimers or for specific interfaces in a multimeric complex
def predict_binding_interface(parent_dir, input_file, pot_dict, z_threshold, output_file, options, p1 = None, p2 = None, f_alascan = None):

	import sys
	import os
	import cPickle as pickle
	from altloc_filter import altloc_check
	from cell_list_main import get_dist
	from parse_dist import parse_dist
	from get_interface_respairs import get_contacts
	from find_clashes_3 import find_clashes
	#from percent_clash import percent_clash, scale_cifa
	from percent_clash_3 import percent_clash, scale_cifa
	from stick_tighter import stick_tighter
	from calc_weight import get_weights
	from get_linkages import get_linkages
	#from score_complex import call_score_func
	from score_complex_clash import call_score_func

	from calc_zscore import calc_zscore
	from alascan import call_alascan
	from get_reslist_als import get_reslist
	from err_warn_report import err_report
	from write_out import write_out

	mutanal_path = os.path.abspath(os.path.join('scripts', 'src_mutanal'))
	sys.path.append(mutanal_path)

	from get_reslist_mutanal import get_mutanal_list
	from mutanal import mutanal_test

	pdb_name = output_file.split('/')[-1].replace('.out', '')

	filtered_pdb = altloc_check(input_file)

	dist_list = get_dist(filtered_pdb, options['cutoff'], p1, p2)

	parsed_dist = parse_dist(dist_list)

	#clashes = find_clashes(parent_dir + '/data/', parsed_dist[0], options['cutoff'])

	warnings = parsed_dist[-1]

	interface_residue_pairs = get_contacts(parsed_dist[0], options['intertype'])

#	inf_res_pairs_all = get_contacts(parsed_dist[0], 'ss')

#	pickle.dump(inf_res_pairs_all[2], open('inf_res_pairs_all.p', 'wb'))

#	n_inf_res_pairs_all = stick_tighter(inf_res_pairs_all[2], 2.0)

#	print n_inf_res_pairs_all

#	for interface in inf_res_pairs_all[2]:
#		for res_pair in inf_res_pairs_all[2][interface]:
#			far_list = inf_res_pairs_all[2][interface][res_pair]
#			near_list = n_inf_res_pairs_all[interface][res_pair]
#			for element in far_list:
#				sp_element = element.split('\t')
#				f_atom_1 = sp_element[0]
#				f_atom_2 = sp_element[1]
#				far_dist = float(sp_element[-1])
#				for item in near_list:
#					if f_atom_1 in item and f_atom_2 in item:
#						near_dist = float(item.split('\t')[-1])
#						if far_dist - near_dist != 2.0:
#							print item, far_dist, near_dist
#							pass

	#clashes = find_clashes(parent_dir + '/data/', inf_res_pairs_all[2], options['cutoff'])

	clashes = find_clashes(parent_dir + '/data/', interface_residue_pairs[2], options['cutoff'])

	#n_clashes = find_clashes(parent_dir + '/data/', n_inf_res_pairs_all, options['cutoff'])

	#for interface in n_clashes:
		#for atom_pair in n_clashes[interface]:
			#print atom_pair, n_clashes[interface][atom_pair]
		#	pass

	clash_file = open(pdb_name + '.clashes', 'w')

	for interface in clashes:
		clash_file.write(interface + ' ' + str(len(clashes[interface])) + '\n')

	scale_dict = percent_clash(interface_residue_pairs[2], clashes)

	num_interactors = interface_residue_pairs[3]

#	num_interactors_all = inf_res_pairs_all[3]

	if len(num_interactors) == 0:
		error_code = 6
		err_report(error_code, output_file)
		sys.exit()

#	if len(num_interactors_all) == 0:
#		error_code = 6
#		err_report(error_code, output_file)
#		sys.exit()

	if options['intertype'] == 'ms':
		weight_dict = get_weights(interface_residue_pairs[2], options, flag_dict = interface_residue_pairs[4])
	else:
		weight_dict = get_weights(interface_residue_pairs[2], options)

#	weight_dict_all = get_weights(inf_res_pairs_all[2], {'intertype' : 'all'})

#--- Ankit was here...
#	s_weight_dict = scale_cifa(weight_dict_all, scale_dict)
	s_weight_dict = scale_cifa(weight_dict, scale_dict)

	scores = call_score_func(parsed_dist[1], interface_residue_pairs[1], s_weight_dict, pot_dict, scale_dict)

#	scores = call_score_func(parsed_dist[1], inf_res_pairs_all[1], s_weight_dict, pot_dict, scale_dict)
#--- Ankit was here...

	native_score = scores[0]
	decoy_scores = scores[1]
	respair_scores = scores[2]

	out_val = calc_zscore(native_score, decoy_scores, output_file)

	linkage_pairs = get_linkages(parsed_dist[0])

	if f_alascan == None:
		pass
	elif len(f_alascan) == 2:
		reslist = get_mutanal_list(f_alascan[1], linkage_pairs, warnings, output_file)
		mutanal = mutanal_test(pot_dict, weight_dict, respair_scores, interface_residue_pairs[0], reslist[0], native_score)
		#reslist = get_reslist(f_alascan[1],linkage_pairs, warnings)
		warnings = reslist[-1]
		#alascan = call_alascan(pot_dict, weight_dict, interface_residue_pairs[0], native_score, mode = f_alascan, res_list = reslist[0])
	else:
		alascan = call_alascan(pot_dict, weight_dict, respair_scores, interface_residue_pairs[0], native_score, mode = f_alascan)

	if f_alascan == None:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], respair_scores, warnings)
	else:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], respair_scores, warnings, alscn = alascan)

	return

