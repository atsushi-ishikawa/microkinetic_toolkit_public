import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
argvs = sys.argv
reactionfile = argvs[1]

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

# r_ads and p_ads are species list of ALL the elementary reactions.
# e.g. if inputfile contains
#  (1) A1 + B1 --> C1 + D1
#  (2) A2 + B2 --> C2 + D2
# it gives 
# r_ads = [ [['A1'],['B1']] , [['A2'],['B2']] ]
# p_ads = [ [['C1'],['D1']] , [['C2'],['D2']] ]

outputfile = "met001ode.m"
fout = open(outputfile,"w")

fout.write('function dydt = ' + outputfile.replace('.m','') + '(~,y,yin,tau,A,Ea,Kci,sden,area,Vr)')
fout.write("\n\n")
#
# template - start
#
template = " \
\t global T R Ngas Ncomp species Ncstr; \n \
\t \n \
\t R  = 8.314; \n \
\t NA = 6.02*10^23; % Avogadro's number \n \
\n \
\t Afor = A(:,1); \n \
\n \
\t kfor = Afor .* exp(-Ea/R/T); \n \
\t krev = kfor ./ Kci; \n \
\n \
\t C = y(1:Ngas);\n  \
\t theta = y(1:Ncomp);\n \
\t theta(1:Ngas) = 0;\n \
\t theta = theta * sden; \n \
\n \
\t Cin = yin(1:Ngas); \n"

fout.write(template)
#
# template - end
#
rxn_num  = get_number_of_reaction(reactionfile)
spec_num = get_species_num()

fout.write("\t Rate = zeros(" + str(spec_num) + ",1);\n\n")

dict1 = {}
dict2 = {}
for irxn in range(rxn_num):
	rxn_idx = str(irxn + 1) # MATLAB-index starts from 1

	# hash-tag based reactant and product species list FOR THIS REACTION
	list_r = []
	list_p = []

	# making dict1, a dictionary with hash of species number and molecule
	for imols, mols in enumerate(r_ads[irxn]):
		for imol, mol in enumerate(mols):
			mol = remove_side_and_flip(mol)
			site = r_site[irxn][imols][imol]
			if site !='gas':
				mol = mol + "_surf"
			spe = get_species_num(mol) + 1 # MATLAB
			list_r.append(spe)
			dict1[spe] = mol

	for imols, mols in enumerate(p_ads[irxn]):
		for imol, mol in enumerate(mols):
			mol = remove_side_and_flip(mol)
			site = p_site[irxn][imols][imol]
			if site !='gas':
				mol = mol + "_surf"
			spe = get_species_num(mol) + 1 # MATLAB
			list_p.append(spe)
			dict1[spe] = mol

	# making dict1 --- done

	#
	# forward direction
	#
	tmp = "kfor(" + rxn_idx + ")"
	for imols, mols in enumerate(r_ads[irxn]):
		for imol, mol in enumerate(mols):
			mol = remove_side_and_flip(mol)
			site = r_site[irxn][imols][imol]

			if site !='gas':
				mol = mol + "_surf"
			spe = get_species_num(mol) + 1 # MATLAB

			if site != 'gas' or mol  == 'surf':
				theta = "theta(" + str(spe) + ")"
			else:
				theta = "C(" + str(spe) + ")"

			power = r_coef[irxn][imols]

			if power != 1:
				theta = theta + "^" + str(power)
			tmp = tmp + "*" + theta

	for mem in list_r:
		coef = 0
		for imol, mol in enumerate(r_ads[irxn]):
			mol = mol[0]
			mol = remove_side_and_flip(mol)

			adsorbate = dict1[mem].split("_")[0]
			if mol == adsorbate:
				coef = r_coef[irxn][imol]
		
		if coef == 0:
 			print("something wrong at coef 1"); exit()

		sto_coef = str(float(coef))
		if mem in dict2:
			dict2[mem] = dict2[mem] + " - " + sto_coef + "*" + tmp
		else:
			dict2[mem] = " - " + sto_coef + "*" + tmp

		# multiply area/Vr for the surface reaction contribution to gas-phase species
		# gas:     k*C = [m^3*s^-1]*[mol*m^-3] = [mol*s^-1]
		# surface: k*(theta*sden)*(A/V)= [m^3*s^-1]*[mol*m^-2]*[m^2/m^3] = [mol*s^-1] --> same!

		if not "surf" in dict1[mem] and "theta" in dict2[mem]:
			dict2[mem] = dict2[mem] + "*" + "(area/Vr)"

	for mem in list_p:
		coef = 0
		for imols, mols in enumerate(p_ads[irxn]):
			for imol, mol in enumerate(mols):
				mol = remove_side_and_flip(mol)
				adsorbate = dict1[mem].split("_")[0]
				if mol == adsorbate:
					coef = p_coef[irxn][imols]
		
		if coef == 0:
			print("something wrong at coef 2"); exit()

		sto_coef = str(float(coef))
		if mem in dict2:
			dict2[mem] = dict2[mem] + " + " + sto_coef + "*" + tmp
		else:
			dict2[mem] = "  " + sto_coef + "*" + tmp

		# multiply area/Vr for the surface reaction contribution to gas-phase species
		if not "surf" in dict1[mem] and "theta" in dict2[mem]:
			dict2[mem] = dict2[mem] + "*" + "(area/Vr)"

	#
	# reverse direction
	#
	tmp = "krev(" + rxn_idx + ")"
	for imols, mols in enumerate(p_ads[irxn]):
		for imol, mol in enumerate(mols):
			mol = remove_side_and_flip(mol)
			site = p_site[irxn][imols][imol]
			if site !='gas':
				mol = mol + "_surf"
			spe = get_species_num(mol) + 1 # MATLAB

			if site != 'gas' or mol  == 'surf':
				theta = "theta(" + str(spe) + ")"
			else:
				theta = "C(" + str(spe) + ")"

			power = p_coef[irxn][imols]

			if power != 1:
				theta = theta + "^" + str(power)
			tmp = tmp + "*" + theta

		# end mol
	# end mols

	for mem in list_r:
		coef = 0
		for imol, mol in enumerate(r_ads[irxn]):
			mol = mol[0]
			mol = remove_side_and_flip(mol)

			adsorbate = dict1[mem].split("_")[0]
			if mol == adsorbate:
				coef = r_coef[irxn][imol]
		
		if coef == 0:
			print("something wrong at coef 3"); exit()

		sto_coef = str(float(coef))
		if mem in dict2:
			dict2[mem] = dict2[mem] + " + " + sto_coef + "*" + tmp
		else:
			dict2[mem] = "  " + sto_coef + "*" + tmp

		# multiply area/Vr for the surface reaction contribution to gas-phase species
		if not "surf" in dict1[mem] and "theta" in dict2[mem]:
			dict2[mem] = dict2[mem] + "*" + "(area/Vr)"

	for mem in list_p:
		coef = 0
		for imols, mols in enumerate(p_ads[irxn]):
			for imol, mol in enumerate(mols):
				mol = remove_side_and_flip(mol)
				adsorbate = dict1[mem].split("_")[0]

				if mol == adsorbate:
					coef = p_coef[irxn][imols]
		
		if coef == 0:
			print("something wrong at coef 4"); exit()

		sto_coef = str(float(coef))
		if mem in dict2:
			dict2[mem] = dict2[mem] + " - " + sto_coef + "*" + tmp
		else:
			dict2[mem] = " - " + sto_coef + "*" + tmp

		# multiply area/Vr for the surface reaction contribution to gas-phase species
		if not "surf" in dict1[mem] and "theta" in dict2[mem]:
			dict2[mem] = dict2[mem] + "*" + "(area/Vr)"
 
#
# end loop for reaction
#

# vacancy site
if 'surf' in dict1.values(): # only when surface is involved
	tmp = ""
	for imol,mol in enumerate(dict1):
		comp = dict1[imol+1]
		if 'surf' in comp and comp != 'surf':
			tmp = tmp + "-Rate(" + str(imol+1) + ")"

	dict2[len(dict2)] = tmp

comment = "\t % species --- "
for imol,mol in enumerate(dict2):
	fout.write("\t Rate({0}) = {1}; % {2}\n".format(imol+1, dict2[imol+1], dict1[imol+1]))
	comment += "%s = %s " % (imol+1,dict1[imol+1])
comment += "\n"

#comment = "\t % species --- " + str(sorted(dict1.items())).replace('\'','').replace('{','').replace('}','') + "\n"
fout.write(comment)
#
# tempelate - start
#
template = "\
\t Rate(Ngas+1:Ncomp) = Rate(Ngas+1:Ncomp)*(area/Vr); % gas + surface \n \
\t % fprintf('------------\\n');\n \
\t % fprintf('%+12.8e %+12.8e %+12.8e %+12.8e %+12.8e \\n', Rate(1),Rate(2),Rate(3),Rate(4),Rate(5));\n \
\t % fprintf('%+12.8e %+12.8e %+12.8e %+12.8e %+12.8e \\n', C(1),C(2),C(3),theta(4),theta(5));\n \
\t if Ncstr == 1 \n \
\t    dydt = Rate(1:Ngas); % gas [mol/m^3]/[sec] = [mol/m^3/sec] \n \
\t else \n \
\t    dydt = 1/tau*(Cin - C) + Rate(1:Ngas); % gas [mol/m^3]/[sec] = [mol/m^3/sec] \n \
\t end \n \
\t if Ncomp > Ngas+1\n \
\t\t dydt(Ngas+1:Ncomp) = Rate(Ngas+1:Ncomp)*(1/sden)*(Vr/area); % species\n \
\t end \n \
"
#
# tempelate - end
#
fout.write(template)

string = ''
#for mol in dict1.values():
for i,_ in enumerate(dict1):
	mol = dict1[i+1]
	if mol == 'surf':
		string = string + '"{0}"'.format('\\theta_{vac}')
	elif 'surf' in mol:
		string = string + '"{0}"'.format('\\theta_{' + mol.replace('_surf','') + '}')
	else:
		string = string + '"{0}"'.format(mol)

	string = string + ","

fout.write( "\t species = [{0}];\n".format(string[:-1]))

fout.write("\nend\n")
fout.close()

