import os
import sys
from reaction_tools import *
#
# make rate equation for python format
#
argvs = sys.argv
reactionfile = argvs[1]

(r_ads, r_site, r_coef, p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

# r_ads and p_ads are species list of ALL the elementary reactions.
# e.g. if inputfile contains
#  (1) A1 + B1 --> C1 + D1
#  (2) A2 + B2 --> C2 + D2
# it gives 
# r_ads = [ [['A1'],['B1']] , [['A2'],['B2']] ]
# p_ads = [ [['C1'],['D1']] , [['C2'],['D2']] ]

outputfile = "ode_equation.py"
fout = open(outputfile, "w")

fout.write('import numpy as np')
fout.write("\n\n")
fout.write('def func(t, c, Afor, Ea, Kci, T, sden, area, Vr, Ngas, Ncomp):')
fout.write("\n")
#
# template - start
#
template = "\
\tR  = 8.314e-3  # kJ/mol/K\n\
\n\
\tkfor = Afor * np.exp(-Ea/R/T)\n\
\tkrev = kfor / Kci\n\
\n\
\ttheta = c[0:Ncomp]\n\
\ttheta = theta * sden\n\
\n"

fout.write(template)
#
# template - end
#
rxn_num = get_number_of_reaction(reactionfile)
spec_num = get_species_num()

fout.write("\trate = np.zeros(" + str(spec_num) + ")\n\n")

dict1 = {}
dict2 = {}
for irxn in range(rxn_num):
	rxn_idx = str(irxn)

	# hash-tag based reactant and product species list FOR THIS REACTION
	list_r = []
	list_p = []
	#
	# making dict1, a dictionary with hash of species-number and molecule
	#
	for side in ["reactant", "product"]:
		if side == "reactant":
			mol_set = r_ads.copy()
			sites   = r_site.copy()
			list    = list_r
		elif side == "product":
			mol_set = p_ads.copy()
			sites   = p_site.copy()
			list    = list_p
		else:
			print("error asdf")
			exit(1)

		for imols, mols in enumerate(mol_set[irxn]):
			for imol, mol in enumerate(mols):
				mol = remove_side_and_flip(mol)
				site = sites[irxn][imols][imol]
				if site != 'gas':
					mol = mol + "_surf"
				spe = get_species_num(mol)
				list.append(spe)
				dict1[spe] = mol
	# done

	for side in ["reactant", "product"]:
		for direction in ["forward", "reverse"]:
			if side == "reactant" and direction == "forward":
				mol_list1 = r_ads.copy()
				sites     = r_site.copy()
				coefs     = r_coef.copy()
				add_to    = list_r
				mol_list2 = r_ads.copy()  # list corresponding to add_to
				term = "kfor[" + rxn_idx + "]"
				sign = " - "
			elif side == "reactant" and direction == "reverse":
				mol_list1 = p_ads.copy()
				sites     = p_site.copy()
				coefs     = p_coef.copy()
				add_to    = list_r
				mol_list2 = r_ads.copy()
				term = "krev[" + rxn_idx + "]"
				sign = " + "
			elif side == "product" and direction == "forward":
				mol_list1 = r_ads.copy()
				sites     = r_site.copy()
				coefs     = r_coef.copy()
				add_to    = list_p
				mol_list2 = p_ads.copy()
				term = "kfor[" + rxn_idx + "]"
				sign = " + "
			elif side == "product" and direction == "reverse":
				mol_list1 = p_ads.copy()
				sites     = p_site.copy()
				coefs     = p_coef.copy()
				add_to    = list_p
				mol_list2 = p_ads.copy()
				term = "krev[" + rxn_idx + "]"
				sign = " - "

			for imols, mols in enumerate(mol_list1[irxn]):
				for imol, mol in enumerate(mols):
					mol  = remove_side_and_flip(mol)
					site = sites[irxn][imols][imol]

					if site != 'gas':
						mol = mol + "_surf"

					spe = get_species_num(mol)
					if site == 'gas':
						if mol == 'surf':  # bare surface
							theta = "theta[" + str(spe) + "]"
						else:  # gas-phase molecule
							theta = "c[" + str(spe) + "]"
					else:  # adsorbed species
						theta = "theta[" + str(spe) + "]"

					power = coefs[irxn][imols]
					if power != 1:
						theta = theta + "**" + str(power)

					term = term + "*" + theta

			for mem in add_to:
				if dict1[mem] == "surf":
					continue

				coef = 0
				for imol, mol in enumerate(mol_list2[irxn]):
					mol = mol[0]
					mol = remove_side_and_flip(mol)

					adsorbate = dict1[mem].split("_")[0]
					if mol == adsorbate:
						coef = coefs[irxn][imol]

				if coef == 0:
					print("something wrong at coef 1")
					exit()

				sto_coef = str(float(coef))

				if mem in dict2:
					dict2[mem] = dict2[mem] + sign + sto_coef + "*" + term  # NEGATIVE
				else:
					dict2[mem] = sign + sto_coef + "*" + term

				if "theta" in dict2[mem]:
					dict2[mem] = dict2[mem] + "*" + "(area/Vr)"


# vacancy site
if 'surf' in dict1.values():  # only when surface is involved
	tmp = ""
	for imol, mol in enumerate(dict1):
		comp = dict1[imol]
		if 'surf' in comp and comp != 'surf':
			tmp = tmp + " -rate[" + str(imol) + "]"

	dict2[len(dict2)] = tmp

comment = "\n\t# species --- "

for imol, mol in enumerate(dict2):
	fout.write("\trate[{0}] ={1}  # {2}\n".format(imol, dict2[imol], dict1[imol]))
	comment += "%s = %s " % (imol, dict1[imol])
comment += "\n"

fout.write(comment)
#
# tempelate - start
#
template = "\
\tif Ncomp > Ngas:\n\
\t\trate[Ngas:Ncomp] = rate[Ngas:Ncomp]*(1/sden)*(Vr/area)  # surface\n\
"
#
# tempelate - end
#
fout.write(template)

# string = ''
# for mol in dict1.values():
# for i, _ in enumerate(dict1):
#	mol = dict1[i]
#	if mol == 'surf':
#		string = string + '"{0}"'.format('\\theta_{vac}')
#	elif 'surf' in mol:
#		string = string + '"{0}"'.format('\\theta_{' + mol.replace('_surf', '') + '}')
#	else:
#		string = string + '"{0}"'.format(mol)
#
#	string = string + ","
#
# fout.write("\t# species = [{0}] \n\n".format(string[:-1]))

fout.write("\treturn rate\n")
fout.close()
