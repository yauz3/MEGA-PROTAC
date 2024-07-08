
def stick_tighter(inf_res_pairs_all, delta_inter_dist):

	n_inf_res_pairs_all = {}

	for interface in inf_res_pairs_all:
		n_inf_res_pairs_all.update({interface : {}})
		for res_pair in inf_res_pairs_all[interface]:
			n_inf_res_pairs_all[interface].update({res_pair : []})
			for atom_pair in inf_res_pairs_all[interface][res_pair]:
				sp_atom_pair = atom_pair.split('\t')
				dist = float(sp_atom_pair[-1])
				n_dist = dist - delta_inter_dist
				if n_dist < 0.0:
					n_dist = 0.0
				n_atom_pair = '\t'.join([sp_atom_pair[0], sp_atom_pair[1], str(n_dist)])
				n_inf_res_pairs_all[interface][res_pair].append(n_atom_pair)

	return n_inf_res_pairs_all
