import os,sys
from reaction_tools import *
#
argvs = sys.argv
reactionfile  = argvs[1]
coveragefile  = argvs[2] # coverage file from MATLAB
rateconstfile = argvs[3] # rate constants from MATLAB
variablefile  = argvs[4] # sden, area, Vr

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
rxn_num  = get_number_of_reaction(reactionfile)
spec_num = get_species_num()

# read coverages
cov = [0 for i in range(spec_num)]
f = open(coveragefile,'r')
lines = f.readlines()
f.close()

for iline, line in enumerate(lines):
	spe = iline
	val = float(line.split()[1])
	cov[spe] = val

# read rate constants
kfor = [0 for i in range(rxn_num)]
krev = [0 for i in range(rxn_num)]
f = open(rateconstfile,'r')
lines = f.readlines()
f.close()

# read variables
f = open(variablefile,'r')
sden = float(f.readline().replace("\n",""))
area = float(f.readline().replace("\n",""))
Vr   = float(f.readline().replace("\n",""))
f.close()

for iline, line in enumerate(lines):
	rxn = iline
	val1,val2 = line.split(' ')
	kfor[iline] = float(val1)
	krev[iline] = float(val2)

rate = [0 for i in range(rxn_num)]

AoverV = area/Vr

for irxn in range(rxn_num):
	#
	# forward
	#
	val = kfor[irxn]
	conc_all = 1.0
	for imol,mol in enumerate(r_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)
		
		site = r_site[irxn][imol]
		if not 'gas' in site:
			# surface species
			mol = mol + '_surf'
			spe = get_species_num(mol)
			conc = cov[spe] * sden * AoverV
		else:
			spe = get_species_num(mol)
			conc = cov[spe]

		power = r_coef[irxn][imol]
		if power != 1:
			conc = conc**power
		conc_all = conc_all * conc

	for_rate = val*conc_all
	#
	# reverse 
	#
	val = krev[irxn]
	conc_all = 1.0
	for imol,mol in enumerate(p_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)

		site = p_site[irxn][imol]
		if not 'gas' in site:
			# surface species
			mol = mol + '_surf'
			spe = get_species_num(mol)
			conc = cov[spe] * sden * AoverV
		else:
			spe = get_species_num(mol)
			conc = cov[spe]

		power = p_coef[irxn][imol]
		if power != 1:
			conc = conc**power
		conc_all = conc_all * conc

	rev_rate = val*conc_all

	rate[irxn] = for_rate - rev_rate

# rate end

reactionfile_out = reactionfile.split('.')[0] + '_val.txt'
f1 = open(reactionfile_out,'w')
lines = return_lines_of_reactionfile(reactionfile)

for iline, line in enumerate(lines):
	text = line.replace('\n','')
	f1.write('{0:<65s}:{1:12.4e}\n'.format(text,rate[iline]))

f1.close()

