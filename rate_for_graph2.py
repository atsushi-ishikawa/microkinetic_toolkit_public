import os, sys
from reaction_tools import *
 
argvs = sys.argv
reactionfile  = argvs[1]
coveragefile  = argvs[2]  # time-dependent coverage file from MATLAB. Saved as "nodes.txt" in MATLAB
rateconstfile = argvs[3]  # rate constants from MATLAB
variablefile  = argvs[4]  # sden, area, Vr

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
rxn_num = get_number_of_reaction(reactionfile)
spe_num = get_species_num()

cov = [[] for i in range(spe_num)]

f = open(coveragefile, 'r')
lines = f.readlines()
lines = lines[1:]
f.close()

Nrxn = 0
for iline, line in enumerate(lines):
	name = str(line.split()[0])
	if "R" not in name:
		spe  = int(line.split()[1])-1  # Matlab --> python
		conc = float(line.split()[2])
		time = float(line.split()[3])

		cov[spe].append([name, conc, time])
	else:
		Nrxn += 1

maxtime = (len(lines)-Nrxn)//spe_num

#
# read rate constants
#
kfor = [0 for i in range(rxn_num)]
krev = [0 for i in range(rxn_num)]

f = open(rateconstfile, 'r')
lines = f.readlines()
f.close()
#
# read variables
#
f = open(variablefile, 'r')
sden = float(f.readline().replace("\n", ""))
area = float(f.readline().replace("\n", ""))
Vr   = float(f.readline().replace("\n", ""))
f.close()

for iline, line in enumerate(lines):
	val1, val2 = line.split(" ")
	kfor[iline] = float(val1)
	krev[iline] = float(val2)

rate = [0 for i in range(rxn_num)]

AoverV = area/Vr

edgefile = 'edges.txt'
f1 = open(edgefile, 'w')
f1.write('    from       to         rate         time\n')

for istep in range(maxtime):
	for irxn in range(rxn_num):
		#
		# reactant side
		#
		val = kfor[irxn]
		conc_all = 1.0
		for imol,mol in enumerate(r_ads[irxn]):
			mol = mol[0]
			mol = remove_side_and_flip(mol)
		
			site = r_site[irxn][imol]
			if 'gas' not in site:
				# surface species
				mol = mol + '_surf'
				spe = get_species_num(mol)
				name, conc, time = cov[spe][istep]
				conc *= sden * AoverV
			else:
				spe = get_species_num(mol)
				name, conc, time = cov[spe][istep]

			power = r_coef[irxn][imol]
			if power != 1:
				conc = conc**power
			conc_all = conc_all * conc

		for_rate = val*conc_all
		#
		# now we write to file
		#
		now = 'R%03d' % irxn # current reaction

		for imol, mol in enumerate(r_ads[irxn]):
			mol = remove_side_and_flip(mol[0])
			site = r_site[irxn][imol]

			if 'gas' not in site:
				mol = mol + '_surf'
			spe = get_species_num(mol)
			name, conc, time = cov[spe][istep]

			f1.write("%8s %8s %12.4e %12.4e\n" % (name, now, for_rate, time))
		
		#
		# product side
		#
		for imol, mol in enumerate(p_ads[irxn]):
			mol = remove_side_and_flip(mol[0])
			site = p_site[irxn][imol]

			if 'gas' not in site:
				mol = mol + '_surf'
			spe = get_species_num(mol)
			name, conc, time = cov[spe][istep]

			f1.write("%8s %8s %12.4e %12.4e\n" % (now, name, for_rate, time))

f1.close()

#f = open(coveragefile, 'a')
#for i in range(rxn_num):
#	f.write("          RXN%03d%5.d%14.5e%14.5f\n" % (i, i+100, 1.0, 0.0))
#f.close()
