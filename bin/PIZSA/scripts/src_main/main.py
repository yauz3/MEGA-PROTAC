"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
Main routine which decides which modules to run based on user input and/or complex type.
'''

def main():
	import os
	import sys
	import time
	try:
		import cPickle as pickle
	except ImportError:
		import pickle
	from parse_pdb import parse_pdb
	from get_options import get_options
	from select_pot import select_pot
	from convert_pot import convert_pot
	from err_warn_report import err_report
	from predict_binding import predict_binding_all, predict_binding_interface

	parent_dir = os.getcwd()
	data_dir = parent_dir + '/data/'
	#input_dir = parent_dir + '/input/'
	#output_dir = parent_dir + '/output/'

	start_time = time.time()

	options = get_options()

	input_dir = '/'.join(options['input_pdb'].split('/')[:-1]) + '/'
	output_dir = options['outfile']

	input_file = options['input_pdb']
	input_file_extn = options['input_pdb'].split('.')[-1]

	if options['outfile'] == None:
		output_file = options['input_pdb'].replace(input_file_extn, 'out')
	else:
		output_file = (options['outfile'] + options['input_pdb'].split('/')[-1]).replace(input_file_extn, '.out')

	parsed_pdb = parse_pdb(input_file)
	chain_set = parsed_pdb[1]
	chain_num = parsed_pdb[0]

	if chain_num < 2:
		error_code = 1
		err_report(error_code, output_file)
		sys.exit()

	operating_points = pickle.load(open(os.path.join(data_dir, 'operating_points.p'), 'rb'))

#--- Ankit was here...
	if options['custom_pot'] == None:
#		pot = select_pot(options['cutoff'], options['intertype'])
#		pot_dict = pickle.load(open(os.path.join(data_dir, pot), 'rb'))
		pot_mcmc = '{}_mcmc_atomprop_matrix.p'.format(int(options['cutoff']))
		pot_mcsc = '{}_mcsc_atomprop_matrix.p'.format(int(options['cutoff']))
		pot_scsc = '{}_scsc_atomprop_matrix.p'.format(int(options['cutoff']))

		pot_dict_mcmc = pickle.load(open(os.path.join(data_dir, pot_mcmc), 'rb'))
		pot_dict_mcsc = pickle.load(open(os.path.join(data_dir, pot_mcsc), 'rb'))
		pot_dict_scsc = pickle.load(open(os.path.join(data_dir, pot_scsc), 'rb'))

		pot_dict = {'mcmc': pot_dict_mcmc,
				'mcsc': pot_dict_mcsc,
				'scsc': pot_dict_scsc}
#--- Ankit was here...

#--- Ankit was here...
		if 'c_ij' not in pot_mcmc:
#			z_threshold = operating_points[pot.replace('.p', '').replace('_', '.')]
			z_threshold = 0.7
		else:
#			z_threshold = operating_points[pot.replace('.p', '').replace('_', '.').replace('c.ij', 'c_ij')]
			z_threshold = - 0.7
#--- Ankit was here...

	else:
		pot_file = open(options['custom_pot'], 'r').readlines()
		pot_dict = convert_pot(pot_file)

		z_threshold = - 0.7

	if options['alascan'] == '0':
		if len(chain_set) > 2:
			if options['protein1'] == None and options['protein2'] == None:
				predict_binding_all(parent_dir, input_file, pot_dict, z_threshold, output_file, options)
			else:
				predict_binding_interface(parent_dir, input_file, pot_dict, z_threshold, output_file, options, p1 = options['protein1'], p2 = options['protein2'])
		else:
#--- Ankit was here...
			predict_binding_all(parent_dir, input_file, pot_dict, z_threshold, output_file, options)
#			predict_binding_interface(parent_dir, input_file, pot_dict, z_threshold, output_file, options, p1 = chain_set[0], p2 = chain_set[1])
#--- Ankit was here...
	else:
		if len(chain_set) > 2:
			if options['protein1'] == None and options['protein2'] == None:
				predict_binding_all(parent_dir, input_file, pot_dict, z_threshold, output_file, options, f_alascan = options['alascan'])
			else:
				predict_binding_interface(parent_dir, input_file, pot_dict, z_threshold, output_file, options, \
				 p1 = options['protein1'], p2 = options['protein2'], f_alascan = options['alascan'])
		else:
#--- Ankit was here...
			predict_binding_all(parent_dir, input_file, pot_dict, z_threshold, output_file, options, f_alascan = options['alascan'])
#			predict_binding_interface(parent_dir, input_file, pot_dict, z_threshold, output_file, options, \
#			 p1 = chain_set[0], p2 = chain_set[1], f_alascan = options['alascan'])
#--- Ankit was here...

	end_time = time.time()

	return
