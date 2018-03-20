import os,sys
from reaction_tools import *
#
# make rate equation for MATLAB format
#
argvs = sys.argv
reactionfile = argvs[1]
#reactionfile = "reaction.txt"

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
outputfile = "met001ode.m"
fout = open(outputfile,"w")

fout.write('function dydt = ' + outputfile.replace('.m','') + '(~,y,yin,tau,A,Ea,deltaH,Wcat,Vr)')
fout.write("\n\n")
#
# template - start
#
template = " \
\t global T R Pin Ngas Ncomp rho_b v0; \n \
\t \n \
\t R = 8.314; \n \
\t NA = 6.02*10^23; % Avogadro's number \n \
\n \
\t Afor = A(:,1); Arev = A(:,2); \n \
\n \
\t kfor = Afor .* exp(-Ea/R/T); \n \
\t krev = Arev .* exp(-(Ea-deltaH)/R/T); \n \
\n \
\t F = y(1:Ngas);\n  \
\t F_tot = sum(F);\n \
\t x = F/F_tot;\n \
\t P = Pin*x;\n \
\t theta = y(1:Ncomp);\n \
\t theta(1:Ngas) = 0;\n \
\n \
\t f_act = 1.0*10^-2; % fraction of active site\n \
\t MW    = 101.1; % atomic mass of Ru [g/mol]\n \
\t Natom = NA/MW; % number of atoms in one g-cat [atoms/g]\n \
\t Nsurf = Natom^(2/3); % number of atoms in surface\n \
\t Nsite = Nsurf*f_act; % number of active sites\n \
\n \
\t Fin = yin(1:Ngas);\n \
\t theta = theta * Nsite;\n \
\t theta = theta * (Wcat/Vr);\n"

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

# vacancy site
if 'surf' in dict1.values(): # only when surface is involved
	tmp = ""
	for imol,mol in enumerate(dict1):
		comp = dict1[imol+1]
		if 'surf' in comp and comp != 'surf':
			tmp = tmp + "-Rate(" + str(imol+1) + ")"

	dict2[len(dict2)] = tmp

for imol,mol in enumerate(dict2):
	fout.write("\t Rate({0}) = {1}; % {2}\n".format(imol+1, dict2[imol+1], dict1[imol+1]))

comment = "\t % species --- " + str(dict1).replace('\'','').replace('{','').replace('}','') + "\n"
fout.write(comment)
#
# tempelate - start
#
template = "\
\t Rate = Rate/NA * Vr;\n \
\t %Rate(1:Ngas)       = Rate(1:Ngas)*(1/NA)*Vr; % [molecule/site/s] --> [mol/s]\n \
\t %Rate(Ngas+1:Ncomp) = Rate(Ngas+1:Ncomp)*(Nsite/NA)*Wcat; % [molecule/site/s] --> [mol/s]\n \
\t %Rate(Ngas+1:Ncomp) = Rate(Ngas+1:Ncomp)*rho_b*v0; % [mol/s] --> [mol/s]*[g/m^3]*[m^3/s]=[mol/s^2]\n \
\n \
\t fprintf('------------\\n');\n \
\t fprintf('%+12.8e %+12.8e %+12.8e %+12.8e %+12.8e \\n', Rate(1),Rate(2),Rate(3),Rate(4),Rate(5));\n \
\t %fprintf('%+12.8e %+12.8e %+12.8e %+12.8e %+12.8e \\n', P(1),P(2),P(3),theta(4),theta(5));\n \
\n \
\t dydt = 1/tau*(F - Fin) + Rate(1:Ngas); % gas\n \
\t dydt(Ngas+1:Ncomp) = Rate(Ngas+1:Ncomp); % species\n \
"
#
# tempelate - end
#
fout.write(template)

string = ''
for mol in dict1.values():
	if mol == 'surf':
		string = string + "'{0}'".format('\\theta_{vac}')
	elif 'surf' in mol:
		string = string + "'{0}'".format('\\theta_{' + mol.replace('_surf','') + '}')
		#string = string + '\theta_{' + mol.replace('_surf','') + '}'
	else:
		string = string + "'{0}'".format(mol)
		#string = string + ' + mol

	string = string + ","

fout.write("\t % legend in plot \n")
fout.write("\t % legend({0})\n".format(string[:-1]))

fout.write("\nend\n")
fout.close()

