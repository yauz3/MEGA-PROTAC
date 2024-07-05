"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine writes the outputs to their respective files.
'''

# Writes out the output files
def write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, respairs, respair_scores, warnings, opt_alpha_fraction, alscn = None):

	import sys

	input_file_extn = options['input_pdb'].split('/')[-1].split('.')[-1]
	pdb_name = options['input_pdb'].replace('pdb','').split('/')[-1].replace(input_file_extn, '')[:-1]

	out_score_file = output_file.replace('.out', '_scores.txt')
	out_respairs_file = output_file.replace('.out', '_respairs.txt')
	out_alascan_file = output_file.replace('.out', '_mut_analysis.txt')
	out_mut_profile_file = output_file.replace('.out', '_mut_profile.txt')
	out_errors_file = output_file.replace('.out', '.err')
	out_warnings_file = output_file.replace('.out', '.warn')

	# Write out the main values (z-score, raw score etc.) in the score file
	sys.stdout = open(out_score_file, 'w')

	for interface in out_val.keys():

#--- Ankit was here...
		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		dist_cutoff = int(options['cutoff'])

		if dist_cutoff == 4:
			z_threshold = 1.72
		elif dist_cutoff == 6:
			z_threshold = 0.92
		else:
			z_threshold = 0.76
#--- Ankit was here...

		print "{0:<40s} {1:1s} {2:<15s}".format('InterfaceBEG', ':', pdb_name + '_' + p1 + '-' + pdb_name + '_' + p2)
		print "{0:<40s} {1:1s} {2:<10s}".format('Input_pdb', ':', pdb_name + '.pdb')
		print "{0:<40s} {1:1s} {2:<15s}".format('Interface', ':', pdb_name + '_' + p1 + '-' + pdb_name + '_' + p2)
		print "{0:<40s} {1:1s} {2:<5d}".format('Number of Interacting Residue Pairs', ':', num_interactors[interface])
		print "{0:<40s} {1:1s} {2:<2.2f}".format('Distance Threshold (A)', ':', options['cutoff'])
		print "{0:<40s} {1:1s} {2:<2s}".format('Potential Type', ':', options['intertype'])
		print "{0:<40s} {1:1s} {2:<6.3f}".format('Raw Score', ':', native_score[interface])
		print "{0:<40s} {1:1s} {2:<6.3f}".format('Normalized Raw Score', ':', native_score[interface] / num_interactors[interface])
		print "{0:<40s} {1:1s} {2:<6.3f}".format('Average Background Score', ':', out_val[interface][1])
		print "{0:<40s} {1:1s} {2:<6.3f}".format('Standard Deviation', ':', out_val[interface][2])
		print "{0:<40s} {1:1s} {2:<3.3f}".format('Z-score', ':', out_val[interface][0])
		print "{0:<40s} {1:1s} {2:<3.3f}".format('Z-score threshold', ':', z_threshold)
		print "{0:<40s} {1:1s} {2:<3.3f}".format('Optimal Alpha fraction', ':', opt_alpha_fraction)
#		print ""

#--- Ankit was here		
#		if out_val[interface][0] < z_threshold:
		if out_val[interface][0] > z_threshold:
#			print "Status: Stable association (It's a BINDER!)"
#			print "{0:<100s}".format("Stable association (BINDER)")
			print "{0:<40s} {1:1s} {2:<100s}".format('Binding Prediction', ':', "Stable association (BINDER)")
		else:
#			print "Status: Unstable association (We have a NON-BINDER!)"
#			print "{0:<100s}".format("Unstable association (NON-BINDER)")
			print "{0:<40s} {1:1s} {2:<100s}".format('Binding Prediction', ':', 'Unstable association (NON-BINDER)')
		print "InterfaceTER"
		print "\n"
#--- Ankit was here



	# Write out the interacting residue pairs for each interface in the respairs file
	sys.stdout = open(out_respairs_file, 'w')

	print "#######Format Description#######\n"
	print "#This file lists all the amino acid residue pairs that constitute the interfaces in the given protein complex."
	print "#The residues are split into columns based on the subunits. Each line represents one interacting pair."
	print "#XXX:YYY:A is the format where XXX = residue number, YYY = residue type & A = protien subunit."
	print "#For example. 321:ALA:A is the Alanine at position 321 of the subunit A."
	print "#######Format Description#######\n"

	for interface in respairs.keys():

		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		print "{0:<10s} {1:1s} {2:<15s}".format('InterfaceBEG', ':', pdb_name + '_' + p1 + '-' + pdb_name + '_' + p2)
		print "\n",

		#print "{0:<15s} {1:4s} {2:<15s} {3:<9s}".format('Residue 1', '   ', 'Residue 2', 'Pair Score')
		#print "\n",

		for res_pair in respairs[interface].keys():
			if res_pair.count('-') == 1:
				res_1 = res_pair.split('-')[0]
				res_2 = res_pair.split('-')[1]
				res_pair = '_'.join([res_1, res_2])
			else:
				split_index = res_pair.find('-')
				if split_index != 0:
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
					res_pair = '_'.join([res_1, res_2])
				else:
					split_index = (res_pair.replace('-', '', 1).find('-')) + 1
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
					res_pair = '_'.join([res_1, res_2])

			pair_score = respair_scores[interface][res_pair]

			resno_1, resno_2 = res_1.split(':')[0], res_2.split(':')[0]
			chain_1, chain_2 = res_1.split(':')[1], res_2.split(':')[1]

			res_pair = res_pair.replace('_', '-')

			resna = respairs[interface][res_pair].split('-')

			resna_1, resna_2 = resna[0], resna[1]

			out_res_1 = resno_1 + ':' + resna_1 + ':' + chain_1
			out_res_2 = resno_2 + ':' + resna_2 + ':' + chain_2

			sort_res_pair = sorted([out_res_1, out_res_2], key = lambda x: x[-1])

			print "{0:<15s} {1:4s} {2:<15s} {3:<6.3f}".format(sort_res_pair[0], "    ", sort_res_pair[1], pair_score)

		print 'InterfaceTER'

		print "\n",


	# Write out the alanine scanning predictions in the alascan file
	if alscn != None:

		if isinstance(alscn.keys()[0], basestring):

			sys.stdout = open(out_alascan_file, 'w')

			print "********************************"
			print "This file contains the results of the mutational analysis."
			print "The most important residues are listed in a ranked order list."
			print "The last column shows the magnitude of change in scores of the protein complexes upon mutagenesis"
#--- Ankit was here...
#			print "Positive values indicate the mutation is destabilising whereas Negative values indicate the mutation is stabilizing"
			print "Positive values indicate the mutation is stabilising whereas Negative values indicate the mutation is destabilizing"
#--- Ankit was here...
			print "********************************"
			print "\n"

			print "{0:<11s} {1:3s} {2:<11s} {3:3s} {4:<11s} {5:3s} {6:<11s}".format('Rank', '   ', 'Residue no.', '   ', 'Chain', '   ', 'Score Diff')

			#print "\t".join(['Rank', 'Residue no.', 'Chain', 'Score Diff'])

			alscn_list = sorted(alscn, key = alscn.get, reverse = True)
			for element in alscn_list:
				rank = alscn_list.index(element) + 1
				res_no = element.split(':')[0]
				chain = element.split(':')[1]
				diff = alscn[element]
				print "{0:<11d} {1:3s} {2:<11s} {3:3s} {4:<11s} {5:3s} {6:<8.3f}".format(rank, '   ', res_no, '   ', chain, '   ', diff)
				#print "\t".join([str(rank), str(res_no), str(chain), str(diff)])

		elif isinstance(alscn.keys()[0], tuple):

			sys.stdout = open(out_alascan_file, 'w')

			print "********************************"
			print "This file contains the results of the mutational analysis."
			print "The first column shows the user-specified residues to be mutated while the second column shows the target mutation."
			print "The last column shows the magnitude of change in scores of the protein complexes upon mutagenesis"
#--- Ankit was here...
#			print "Positive values indicate the mutation is destabilising whereas Negative values indicate the mutation is stabilising"
			print "Positive values indicate the mutation is stabilising whereas Negative values indicate the mutation is destabilising"
#--- Ankit was here...
			print "********************************"
			print "\n"

			print "{0:<15s} {1:3s} {2:<15s} {3:3s} {4:<15s} {5:3s} {6:<15s} {7:3s} {8:<15s}".format('Residue number', '   ', 'Residue chain', '   ', 'Residue type', '   ', 'Mutant Type', '   ', 'Score Diff')

			mutate_to_set = set()
			for element in alscn:
				mutate_to_set.add(element[1])

			if len(mutate_to_set) > 1:
				for element in alscn:
					split_res_in = element[0].split(':')
 					print "{0:<15s} {1:3s} {2:<15s} {3:3s} {4:<15s} {5:3s} {6:<15s} {7:3s} {8:<12.3f}".format(split_res_in[0], '   ', split_res_in[2], '   ', split_res_in[1], '   ', element[1], '   ', alscn[element][0])
 			else:
 				outlines = []
 				for element in alscn:
 					split_res_in = element[0].split(':')
 					out_string = [split_res_in[0], '   ', split_res_in[2], '   ', split_res_in[1], '   ', element[1], '   ', alscn[element][0]]
 					outlines.append(out_string)

 				outlines = sorted(outlines, key = lambda x: x[-1], reverse = True)

 				for line in outlines:
 					print "{0:<15s} {1:3s} {2:<15s} {3:3s} {4:<15s} {5:3s} {6:<15s} {7:3s} {8:<12.3f}".format(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8])


 			sys.stdout = open(out_mut_profile_file, 'w')

 			print "********************************"
 			print "This file contains the results of the mutational analysis."
 			print "The first column shows the user-specified residues to be mutated while the second column shows the target mutation."
 			print "The third column represents the interface being impacted due to the mutation."
 			print "The last column shows the magnitude of change in scores of the protein complexes upon mutagenesis"
#--- Ankit was here...
# 			print "Positive values indicate the mutation is destabilsing whereas Negative values indicate the mutation is stabilizing"
 			print "Positive values indicate the mutation is stabilsing whereas Negative values indicate the mutation is destabilizing"
#--- Ankit was here...
 			print "********************************"
 			print "\n"

 			print "{0:<11s} {1:3s} {2:<11s} {3:3s} {4:<11s} {5:3s} {6:<11s}".format('Input Res', '   ', 'Mutant Type', '   ', 'Interface', '   ', 'Score Diff')

 			for element in alscn:
 				if element[0].split(':')[1] != element[1]:
 					for interface in alscn[element][1]:
 						print "{0:<11s} {1:3s} {2:<11s} {3:3s} {4:<11s} {5:3s} {6:<8.3f}".format(element[0], '   ', element[1], '   ', interface, '   ', alscn[element][1][interface])



	# Write out the warnings (if any)
	if len(warnings) != 0:

		sys.stdout = open(out_warnings_file, 'w')

		from collections import OrderedDict

		warn_codes_dict = OrderedDict()

		for element in warnings:
			warn_code = element[0]
			if warn_code not in warn_codes_dict:
				warn_codes_dict.update({warn_code : [element[1]]})
			else:
				warn_codes_dict[warn_code].append(element[1])

		for code in warn_codes_dict:

			res_seen = []

			if warn_code == 1:
				print "Code 1: Nucleotides found on the interface \n"
				for res in warn_codes_dict[warn_code]:
					res = ':'.join(res.split(':')[1:])
					if res not in res_seen:
						print res
						res_seen.append(res)

			elif warn_code == 2:
				print "Code 2: Non-Standard Amino Acids found on the interface \n"
				for res in warn_codes_dict[warn_code]:
					res = ':'.join(res.split(':')[1:])
					if res not in res_seen:
						print res
						res_seen.append(res)
			
			elif warn_code == 3:
				print "Code 3: Non-interface residues provided for mutational analysis \n"
				for res in warn_codes_dict[warn_code]:
					res = ':'.join(res.split(':')[1:])
					if res not in res_seen:
						print res
						res_seen.append(res)



	# Write out the errors (if any) (This is implemented in the err_warn_report Module)
	#if error_code != 0:

	#	sys.stdout = open(out_errors_file, 'w')

	#	if error_code == 1:
	#		print "Error Code 1: \nInput file does not contain a Protein Complex"

	#	elif error_code == 2:
	#		print "Error Code 2: \nStandard Deviation of background set is 0"

	#	elif error_code == 3:
	#		print "Error Code 3: \nPlease provide a Residue list for Mutational Analysis"




	# Write out the respair files for our webserver developer
	for interface in respairs.keys():
	
		sys.stdout = open(output_file.replace('.out', '_' + interface + '_interface_respairs.info'), 'w')

		interface_res = set()

		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		for res_pair in respairs[interface].keys():
			if res_pair.count('-') == 1:
				res_1 = res_pair.split('-')[0]
				res_2 = res_pair.split('-')[1]
				res_pair = '_'.join([res_1, res_2])
			else:
				split_index = res_pair.find('-')
				if split_index != 0:
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
					res_pair = '_'.join([res_1, res_2])
				else:
					split_index = (res_pair.replace('-', '', 1).find('-')) + 1
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
					res_pair = '_'.join([res_1, res_2])

			pair_score = respair_scores[interface][res_pair]

			resno_1, resno_2 = res_1.split(':')[0], res_2.split(':')[0]
			chain_1, chain_2 = res_1.split(':')[1], res_2.split(':')[1]

			res_pair = res_pair.replace('_', '-')

			resna = respairs[interface][res_pair].split('-')

			resna_1, resna_2 = resna[0], resna[1]

			out_res_1 = resno_1 + ':' + resna_1 + ':' + chain_1
			out_res_2 = resno_2 + ':' + resna_2 + ':' + chain_2

			interface_res.add(out_res_1)
			interface_res.add(out_res_2)

		for element in interface_res:
			print "{0:<15s}".format(element)


		sys.stdout = open("/dev/stdout", "w")

	# Per residue scores for Neeladri
	#sys.stdout = open(output_file.replace('.out', '_res_scores.out'), 'w')

	#for residue in res_scores.keys():
	#	print "{0:<10s} {1:10.3f}".format(residue, res_scores[residue])

 
	#sys.stdout = open("/dev/stdout", "w")

	return
