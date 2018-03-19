import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
argvs = sys.argv
reactionfile = argvs[1]
#reactionfile = "reaction.txt"

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
outputfile = "met001ode_auto.m"
fout = open(outputfile,"w")

# Remove "surf" from the list, because ODE for vacant site is 
# determined from site balance i.e. sum_i theta_i = 1.
# Note that not necessary remove from species list.
for lst in [r_ads, p_ads]:
	for ads in lst:
		if 'surf' in ads:
			ads.remove('surf')

fout.write("function dydt = met001ode_auto(~,y,yin,tau,A,Ea,deltaH)")
fout.write("\n\n")
#
# template zone - start
#
template = " \
\t global T R Pin Ngas rho_b v0; \n \
\t \n \
\t R = 8.314; \n \
\t NA = 6.02*10^23; % Avogadro's number \n \
\t \n \
\t Afor = A(:,1); Arev = A(:,2); \n \
\t \n \
\t kfor = Afor .* exp(-Ea/R/T); \n \
\t krev = Arev .* exp(-(Ea-deltaH)/R/T); \n \
\t \n \
\t F = y(1:Ngas); \n  \
\t F_tot = sum(F); \n \
\t x = F/F_tot; \n \
\t P = Pin*x; \n \
\t theta = y(Ngas+1:end); \n \
\t \n \
\t f_act = 1.0e-1; % fraction of active site \n \
\t MW    = 101.1; % atomic mass of Ru [g/mol] \n \
\t Natom = NA/MW; % number of atoms in one g-cat [atoms/g] \n \
\t Nsurf = Natom^(2/3); % number of atoms in surface \n \
\t Nsite = Nsurf*f_act; % number of active sites \n \
\t \n \
\t Fin = yin(1:Ngas);\n"
fout.write(template)
#
# template zone - end
#
rxn_num  = get_number_of_reaction(reactionfile)
spec_num = get_species_num()

fout.write("\t Rate = zeros(" + str(spec_num+1) + ",1);\n\n")  # add +1 for vacant site

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
	tmp = "kfor(" + rxn_idx + ")"
	for imol,mol in enumerate(r_ads[irxn]):
		site = r_site[irxn][imol]
		if site !='gas':
			mol = mol + "_surf"
		spe  = get_species_num(mol) + 1 # MATLAB

		if site != 'gas' or mol  == 'surf':
			theta = "theta(" + str(spe) + ")"
		else:
			theta = "P(" + str(spe) + ")"

		power = r_coef[irxn][imol]
		if power != 1:
			theta = theta + "^" + str(power)
		tmp = tmp + "*" + theta

	for mem in list_r:
		if mem in dict2:
			dict2[mem] = dict2[mem] + "-" + tmp
		else:
			dict2[mem] = "-" + tmp
	for mem in list_p:
		if mem in dict2:
			dict2[mem] = dict2[mem] + "+" + tmp
		else:
			dict2[mem] = tmp

	# backword reaction
	tmp = "krev(" + rxn_idx + ")"
	for imol,mol in enumerate(p_ads[irxn]):
		site = p_site[irxn][imol]
		if site !='gas':
			mol = mol + "_surf"
		spe = get_species_num(mol) + 1 # MATLAB

		if site != 'gas' or mol  == 'surf':
			theta = "theta(" + str(spe) + ")"
		else:
			theta = "P(" + str(spe) + ")"

		power = p_coef[irxn][imol]
		if power != 1:
			theta = theta + "^" + str(power)
		tmp = tmp + "*" + theta

	for mem in list_r:
		if mem in dict2:
			dict2[mem] = dict2[mem] + "+" + tmp
		else:
			dict2[mem] = tmp
	for mem in list_p:
		if mem in dict2:
			dict2[mem] = dict2[mem] + "-" + tmp
		else:
			dict2[mem] = "-" + tmp


for imol,mol in enumerate(dict2):
	fout.write("\t Rate({0}) = {1}; % {2}\n".format(imol+1, dict2[imol+1], dict1[imol+1]))

comment = "\t % species --- " + str(dict1).replace('\'','').replace('{','').replace('}','')
fout.write(comment)

fout.write("\nend\n")
fout.close()

