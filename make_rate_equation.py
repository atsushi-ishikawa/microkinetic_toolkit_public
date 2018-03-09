import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
argvs = sys.argv
reactionfile = argvs[1]
#reactionfile = "reaction.txt"

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
outputfile = "dtheta_dt.m"
fout = open(outputfile,"w")

# remove "surf" from the list
for i in r_ads:
	try:
		i = i.remove('surf')
	except:
		i = i
for i in p_ads:
	try:
		i = i.remove('surf')
	except:
		i = i

fout.write("function dtheta_dt = f(t,theta,kfor,krev,P)")
fout.write("\n\n")
#
# template zone - start
#

#
# template zone - end
#
rxn_num = get_number_of_reaction(reactionfile)

fout.write("Rate = zeros(" + str(rxn_num) + ",1);\n\n")

dict = {}
for irxn in range(rxn_num):
	rxn_idx = str(irxn + 1) # MATLAB-index starts from 1

	tmp = "Rate(" + rxn_idx + ") = "
	#
	# reactants
	#
	tmp = tmp + "kfor(" + rxn_idx + ")"

	for imol, mol in enumerate(r_ads[irxn]):
		spe = get_species_num(mol)
		spe_idx = str(spe + 1)  # MATLAB-index
		dict[spe+1] = mol

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
		rxn_idx = str(irxn + 1) # MATLAB-index
		spe = get_species_num(mol)
		spe_idx = str(spe + 1)  # MATLAB-index
		dict[spe+1] = mol

		if p_site[irxn][imol] == "gas":
			theta = "P(" + spe_idx + ")"
		else:
			theta = "theta(" + spe_idx + ")"

		power = p_coef[irxn][imol]
		if power != 1:
			theta = theta + "^" + str(power)
		tmp = tmp + "*" + theta

	print "for and rev", tmp
	tmp = tmp + "\n"
	fout.writelines(tmp)

comment = "% " + str(dict).replace('\'','').replace('{','').replace('}','')
fout.write(comment)

fout.write("\nend\n")
fout.close()


