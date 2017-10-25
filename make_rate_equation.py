import os,sys
from reaction_tools import *
#
# calculate reaction energy
# molecule's data should be stored in "methane.json"
#
reactionfile = "reaction.txt"
(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

rxn_num = get_number_of_reaction(reactionfile)

for irxn in range(rxn_num):
	rxn_idx = str(irxn + 1) # MATLAB

	tmp = "dtheta_dt(" + rxn_idx + ") = "
	#
	# reactants
	#
	tmp = tmp + "kfor(" + rxn_idx + ")"

	for imol, mol in enumerate(r_ads[irxn]):
		spe = get_species_num(mol)
		spe_idx = str(spe + 1)  # MATLAB

		if r_site[irxn][imol] == "gas":
			theta = "P(" + rxn_idx + ")"
		else:
			theta = "theta(" + spe_idx + ")"
		tmp = tmp + "*" + theta

	#
	# products
	#
	tmp = tmp + " - krev(" + rxn_idx + ")"

	for imol, mol in enumerate(p_ads[irxn]):
		rxn_idx = str(irxn + 1) # MATLAB
		spe = get_species_num(mol)
		spe_idx = str(spe + 1)  # MATLAB

		if p_site[irxn][imol] == "gas":
			theta = "P(" + rxn_idx + ")"
		else:
			theta = "theta(" + spe_idx + ")"
		tmp = tmp + "*" + theta


	print "for and rev", tmp

