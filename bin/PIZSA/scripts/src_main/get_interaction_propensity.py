
#--- Interface propensity calculation module of PIZSA alpha 1.7

import os
import cPickle as pickle
import numpy as np
from scipy.interpolate import interp1d


def getPropensity(distances, dist_cutoff, intertype):


#-------#--- Choose allowed atoms
	def chooseAtoms(intertype):

		# allowed side-chain atoms
		goodAtoms = ['CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2', 'CZ', 'CZ2', 'CZ3',
				'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
				'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH',
				'SD', 'SG']

		# add main chain atoms when intertype is all
		if intertype == 'all':
			goodAtoms.extend(['C', 'CA', 'N', 'O'])

		return goodAtoms


#-------#--- Get atoms in every interface residue pair
	def getIntRes(distances):
		
		mcatoms = ['C', 'CA', 'N', 'O']		# main chain atoms

		interface_residues = {}			# stores the atoms involved in each residue pair interaction

		resID_track = set()			# keeps a track of previously encountered residue ID

		#--- Allowed atoms
		goodAtoms = chooseAtoms(intertype)	# get allowed atoms

		#--- Get atoms involved in residue pair interaction from all interfaces
		for interface in distances.keys():
			interactions = distances[interface]		# all interactions in an interface

			for line in interactions:
				part1 = line.split('\t')[0]		# cell list output column 1
				part2 = line.split('\t')[1]		# cell list output column 2

				atom1 = part1.split(':')[0]		# atom from residue 1
				atom2 = part2.split(':')[0]		# atom from residue 2


				#--- Keep interaction pair only if allowed atoms participate in the interaction
				if (atom1 in goodAtoms) and (atom2 in goodAtoms):
					residP1 = ':'.join(part1.split(':')[1:])
					residP2 = ':'.join(part2.split(':')[1:])

					#--- Create resID:
					#    Res1:ResNum:Chain-Res2:ResNum:Chain
					resID = '-'.join([residP1, residP2])

				else:
					continue


				#--- If residue ID is encountered for the first time, create a new entry
				#    Else add to pre-existing residue ID
				if resID not in resID_track:

					#--- Create new entry with MC-MC, MC-SC and SC-SC fields
					interface_residues[resID] = {'mcmc': set(), 'mcsc': set(), 'scsc': set()}

					resID_track.add(resID)		# track resID

					if (atom1 in mcatoms) and (atom2 in mcatoms):
						interface_residues[resID]['mcmc'].add(part1)	# add mcmc interactor 1
						interface_residues[resID]['mcmc'].add(part2)	# add mcmc interactor 2

					elif (atom1 in mcatoms) or (atom2 in mcatoms):
						interface_residues[resID]['mcsc'].add(part1)	# add mcsc interactor 1
						interface_residues[resID]['mcsc'].add(part2)	# add mcsc interactor 2

					else:
						interface_residues[resID]['scsc'].add(part1)	# add scsc interactor 1
						interface_residues[resID]['scsc'].add(part2)	# add scsc interactor 2

				else:

					if (atom1 in mcatoms) and (atom2 in mcatoms):
						interface_residues[resID]['mcmc'].add(part1)	# add mcmc interactor 1
						interface_residues[resID]['mcmc'].add(part2)	# add mcmc interactor 2

					elif (atom1 in mcatoms) or (atom2 in mcatoms):
						interface_residues[resID]['mcsc'].add(part1)	# add mcsc interactor 1
						interface_residues[resID]['mcsc'].add(part2)	# add mcsc interactor 2

					else:
						interface_residues[resID]['scsc'].add(part1)	# add scsc interactor 1
						interface_residues[resID]['scsc'].add(part2)	# add scsc interactor 2


		combined_IntRes = {}
		combined_IntRes['X-X'] = interface_residues

		return combined_IntRes
		


#-------#--- Count atoms participating in a residue pair interaction
	def countAtoms(interface_residues):

		atomCount = {}		# stores the count of atoms in every residue pair interaction

		#--- Get atom count for all residue pair interactions in all interfaces
		for interface in interface_residues.keys():
			count = {}	# stores atom count for all resID

			for resID in interface_residues[interface]:
				count[resID] = {}

				count[resID]['mcmc'] = len(interface_residues[interface][resID]['mcmc'])
				count[resID]['mcsc'] = len(interface_residues[interface][resID]['mcsc'])
				count[resID]['scsc'] = len(interface_residues[interface][resID]['scsc'])

			atomCount[interface] = count

		return atomCount



#-------#--- Load Normalized Propensities
	def loadAtomCountProp():

		parentDir = os.getcwd()
		dataDir = parentDir + '/data/'				# data directory

		atomcountprop_mcmc_filename = '{}_mcmc_atomcount_propensity.p'.format(int(dist_cutoff))
		atomcountprop_mcsc_filename = '{}_mcsc_atomcount_propensity.p'.format(int(dist_cutoff))
		atomcountprop_scsc_filename = '{}_scsc_atomcount_propensity.p'.format(int(dist_cutoff))

		#--- Load normalized propensities
		atomCountProp_mcmc = pickle.load(open(os.path.join(dataDir, atomcountprop_mcmc_filename), 'rb'))	# mcmc propensity
		atomCountProp_mcsc = pickle.load(open(os.path.join(dataDir, atomcountprop_mcsc_filename), 'rb'))	# mcsc propensity
		atomCountProp_scsc = pickle.load(open(os.path.join(dataDir, atomcountprop_scsc_filename), 'rb'))	# scsc propensity

		#--- All normalized propensities
		atomCountProp = {'mcmc': atomCountProp_mcmc,
				 'mcsc': atomCountProp_mcsc,
				 'scsc': atomCountProp_scsc}

		return atomCountProp



#-------#--- Interpolate propensity
#	def interpolation(normProp, obs_respair, obs_count):

#		x = np.array(normProp[obs_respair].keys())		# atom counts
#		y = np.array(normProp[obs_respair].values())		# propensity of atom counts

#		minx = np.min(x)					# minimum atoms
#		maxx = np.max(x)					# maximum atoms

#		lessdata = False					# turns True when less than or equal to 2 data points encountered for interpolation

#		try:
#			ifunc = interp1d(x, y, kind='quadratic', bounds_error=False)		# fitting function
#		except:
#			lessdata = True					# insufficient data to interpolate

		#--- Find propensity
		#--- Interpolate when enough data points and observed count within bounds
#		if ((obs_count >= minx) and (obs_count <= maxx)) and not lessdata:
#			propensity = float(ifunc(obs_count))

		#--- Get raw value when less data but within bounds
#		elif (obs_count >= minx) and (obs_count <= maxx):
#			propensity = normProp[obs_respair][round(obs_count)]

		#--- Zero when out of bounds
#		else:
#			propensity = 0


#		return propensity



#-------#--- Split residue IDs; handles both positive and negative integers as residue numbers
	def resid_splitter(resID):

		negres = resID.count('-') > 1

		#--- At least one negative residue number
		if negres:
			parts = resID.split(':')		# Res, ResNum, Chain-Res, ResNum, Chain
			center = parts[2]			# Chain-Res
			
			resid1 = ':'.join([parts[1], center.split('-')[0]])
			resid2 = ':'.join([parts[3], parts[4]])

		#--- All positive residue numbers
		else:
			resid1 = ':'.join(resID.split('-')[0].split(':')[1:])
			resid2 = ':'.join(resID.split('-')[1].split(':')[1:])

		return (resid1, resid2)



#-------#--- Calculate propensities
	def calcProp(atomCount, atomCountProp):

		interface_prop = {}				# store interaction propensities of every residue pair in all interfaces

		tot_int_count = 0				# total interaction count
		tot_opt_count = 0				# total number of optimal interactions

		for interface in atomCount.keys():
			interface_prop[interface] = {}		# store interaction propensity of residue pairs in an interface

			for resID in atomCount[interface].keys():
				res1 = resID.split('-')[0].split(':')[0]	# residue 1
				res2 = resID.split('-')[1].split(':')[0]	# residue 2

				obs_respair = '-'.join([res1, res2])		# residue pair

				obs_count_mcmc = atomCount[interface][resID]['mcmc']		# number of atoms involved in mcmc interaction
				obs_count_mcsc = atomCount[interface][resID]['mcsc']		# number of atoms involved in mcsc interaction
				obs_count_scsc = atomCount[interface][resID]['scsc']		# number of atoms involved in scsc interaction

				try:
					propensity_mcmc = atomCountProp['mcmc'][obs_respair][obs_count_mcmc]		# propensity of mcmc interaction

					#--- Increment on encountering interaction
					if propensity_mcmc != 'NA':
						tot_int_count += 1

					#--- Increment on encountering optimal interaction
					if propensity_mcmc > 1 and propensity_mcmc != 'NA':
						tot_opt_count += 1
				except:
					propensity_mcmc = 0								# set propensity to 0 for illegal counts

				try:
					propensity_mcsc = atomCountProp['mcsc'][obs_respair][obs_count_mcsc]		# propensity of mcsc interaction

					#--- Increment on encountering interaction
					if propensity_mcsc != 'NA':
						tot_int_count += 1

					#--- Increment on encountering optimal interaction
					if propensity_mcsc > 1 and propensity_mcsc != 'NA':
						tot_opt_count += 1
				except:
					propensity_mcsc = 0								# set propensity to 0 for illegal counts

				try:
					propensity_scsc = atomCountProp['scsc'][obs_respair][obs_count_scsc]		# propensity of scsc interaction

					#--- Increment on encountering interaction
					if propensity_scsc != 'NA':
						tot_int_count += 1

					#--- Increment on encountering optimal interaction
					if propensity_scsc > 1 and propensity_scsc != 'NA':
						tot_opt_count += 1
				except:
					propensity_scsc = 0								# set propensity to 0 for illegal counts

				resid1, resid2 = resid_splitter(resID)
#				resid1 = ':'.join(resID.split('-')[0].split(':')[1:])
#				resid2 = ':'.join(resID.split('-')[1].split(':')[1:])
				newResID = '_'.join([resid1, resid2])				# new residue ID --> ResNum:Chain_ResNum:Chain

				propensity = {'mcmc': propensity_mcmc,
						'mcsc': propensity_mcsc,
						'scsc': propensity_scsc}

				interface_prop[interface][newResID] = propensity		# store propensity

		#--- Fraction of interactions with optimal alphas
		try:
			opt_alpha_fraction = float(tot_opt_count) / float(tot_int_count)
		except:
			opt_alpha_fraction = 0

		return interface_prop, opt_alpha_fraction



#-------#--- Control flow
	def control():
		combined_IntRes = getIntRes(distances)					# get unique interacting atoms in every residue pair interaction
		atomCount = countAtoms(combined_IntRes)					# count the number of interacting atoms
		atomCountProp = loadAtomCountProp()					# load normalized propensities for every residue pair
		interface_prop, opt_alpha_fraction = calcProp(atomCount, atomCountProp)	# get propensities for all residue pair interactions

		return interface_prop, opt_alpha_fraction


	interface_prop, opt_alpha_fraction = control()

	return interface_prop, opt_alpha_fraction
