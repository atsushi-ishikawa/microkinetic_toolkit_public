import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
reactionfile = "reaction.txt"
(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
outputfile = "dtheta_dt.m"
fout = open(outputfile,"w")

fout.write("function dtheta_dt = f(t,theta,kfor,krev,P)")
fout.write("\n\n")

rxn_num = get_number_of_reaction(reactionfile)

fout.write("dtheta_dt = zeros(" + str(rxn_num) + ",1);\n\n")

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
			theta = "P(" + spe_idx + ")"
		else:
			theta = "theta(" + spe_idx + ")"

		power = r_coef[irxn][imol]
		if power != 1:
			theta = theta + "^" + str(power)

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
			theta = "P(" + spe_idx + ")"
		else:
			theta = "theta(" + spe_idx + ")"

		power = p_coef[irxn][imol]
		if power != 1:
			theta = theta + "^" + str(power)
		tmp = tmp + "*" + theta

	tmp = tmp + "\n"
	print "for and rev", tmp
	fout.writelines(tmp)

fout.write("\nend\n")
fout.close()

