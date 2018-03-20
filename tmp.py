import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
argvs = sys.argv
reactionfile = argvs[1]

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
rxn_num  = get_number_of_reaction(reactionfile)
spec_num = get_species_num()

dict1 = {}
dict2 = {}
for irxn in range(rxn_num):
	rxn_idx = str(irxn + 1) # MATLAB-index starts from 1

	list_r = []
	list_p = []

	# prepare dict1 -- start
	for imol,mol in enumerate(r_ads[irxn]):
		site = r_site[irxn][imol]
		if site !='gas':
			mol = mol + "_surf"
		spe  = get_species_num(mol) + 1 # MATLAB
		list_r.append(spe)
		dict1[spe] = mol

	for imol,mol in enumerate(p_ads[irxn]):
		site = p_site[irxn][imol]
		if site !='gas':
			mol = mol + "_surf"
		spe  = get_species_num(mol) + 1 # MATLAB
		list_p.append(spe)
		dict1[spe] = mol
	# prepare dict1 -- end

	#
	# forward reaction
	#
	val = kfor[rxn_idx]
	for imol,mol in enumerate(r_ads[irxn]):
		spe   = get_species_num(mol)
		conc  = cov[spe]
		power = r_coef[irxn][imol]
		if power != 1:
			conc = conc**power
