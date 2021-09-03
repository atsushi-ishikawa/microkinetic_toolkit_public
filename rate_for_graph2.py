import numpy as np
import os, sys
import argparse
import pickle
from reaction_tools import *

parser = argparse.ArgumentParser()
parser.add_argument("--reactionfile", required=True, help="file with elementary reactions")
parser.add_argument("--coveragefile", default="nodes.txt", help="time-dependent coverage")
argvs = parser.parse_args()

reactionfile  = argvs.reactionfile
coveragefile  = argvs.coveragefile
rateconstfile = "rateconst.pickle"
variablefile  = "variables.pickle"   # sden, area, Vr

(r_ads, r_site, r_coef, p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)
rxn_num = get_number_of_reaction(reactionfile)
spe_num = get_species_num()

Nrxn = 0
#
# read coverage
#
cov = [[] for i in range(spe_num)]

with open(coveragefile, 'r') as f:
	lines = f.readlines()
	lines = lines[1:]

for iline, line in enumerate(lines):
	name = str(line.split()[0])
	if "R" not in name:
		spe = int(line.split()[1]) - 1  # Matlab --> python
		conc = float(line.split()[2])
		time = float(line.split()[3])

		cov[spe].append([name, conc, time])
	else:
		Nrxn += 1

max_time = (len(lines) - Nrxn) // spe_num

# read rate constants
kfor, krev = pickle.load(open(rateconstfile, "rb"))

# read variables
sden, area, Vr = pickle.load(open(variablefile, "rb"))

AoverV = area / Vr
#
# output information for graph
#
rates = {"forward": None, "reverse": None}
for direction in ["forward", "reverse"]:
	if direction == "forward":
		mols  = r_ads.copy()
		sites = r_site.copy()
		coefs = r_coef.copy()
		k = kfor.copy()
	elif direction == "reverse":
		mols  = p_ads.copy()
		sites = p_site.copy()
		coefs = p_coef.copy()
		k = krev.copy()

	rate = np.zeros((max_time, rxn_num))
	for istep in range(max_time):
		for irxn in range(rxn_num):

			conc_all = 1.0
			for imol, mol in enumerate(mols[irxn]):
				mol = mol[0]
				mol = remove_side_and_flip(mol)
				site = sites[irxn][imol]

				if 'gas' not in site:  # surface species
					mol += '_surf'
					spe = get_species_num(mol)
					name, conc, time = cov[spe][istep]
					conc *= sden * AoverV
				else:  # gas
					spe = get_species_num(mol)
					name, conc, time = cov[spe][istep]

				power = coefs[irxn][imol]

				if power != 1:
					conc = conc ** power

				conc_all = conc_all * conc

			rate[istep][irxn] = k[irxn] * conc_all

		rates[direction] = rate

#
# now write to file
#
# time-dependent case
#
edgefile = "edges.txt"
f = open(edgefile, "w")
f.write("    from       to         rate         time\n")

for direction in ["forward", "reverse"]:
	if direction == "forward":
		mols = r_ads.copy()
		sites = r_site.copy()
		coefs = r_coef.copy()
	elif direction == "reverse":
		mols = p_ads.copy()
		sites = p_site.copy()
		coefs = p_coef.copy()

	for istep in range(max_time):
		for irxn in range(rxn_num):
			now = 'R%03d' % irxn  # current reaction

			for imol, mol in enumerate(mols[irxn]):
				mol = remove_side_and_flip(mol[0])
				site = sites[irxn][imol]

				if 'gas' not in site:
					mol = mol + '_surf'

				spe = get_species_num(mol)
				name, conc, time = cov[spe][istep]

				if direction == "forward":
					f.write("%10s %10s %12.4e %12.4e\n" % (name, now, rates[direction][istep][irxn], time))
				elif direction == "product":
					f.write("%10s %10s %12.4e %12.4e\n" % (now, name, rates[direction][istep][irxn], time))

f.close()

# time-independent case

reactionfile_out = reactionfile.split('.')[0] + '_val.txt'
f = open(reactionfile_out, 'w')
lines = return_lines_of_reactionfile(reactionfile)

for iline, line in enumerate(lines):
	text = line.replace('\n', '')
	total_rate = rates["forward"][-1][iline] - rates["reverse"][-1][iline]
	f.write('{0:<80s}:{1:12.4e}\n'.format(text, total_rate))

f.close()
