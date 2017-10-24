import os,sys
from reaction_tools import *
#
# calculate reaction energy
# molecule's data should be stored in "methane.json"
#
reactionfile = "reaction.txt"
(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

rxn_num = get_number_of_reaction(reactionfile)

tmp = "delta_dt(1) = "
for irxn in range(rxn_num):
	print "---", irxn, "---"
	#
	# reactants
	#
	tmp = tmp + "kfor[" + irxn + "]"
	for imol, mol in enumerate(r_ads[irxn]):
		num = get_species_num(mol)
		if r_site[irxn][imol] == "gas":
			theta = "P[" + num + "]"
		else
			theta = "theta[" + num + "]"

		print "forward", tmp

	#
	# products
	#
	tmp = tmp + "- krev[" + irxn + "]"
	for imol, mol in enumerate(p_ads[irxn]):
		num = get_species_num(mol)
		if p_site[irxn][imol] == "gas":
			theta = "P[" + num + "]"
		else
			theta = "theta[" + num + "]"

	print "for and rev",tmp

