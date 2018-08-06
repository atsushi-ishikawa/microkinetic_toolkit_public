import numpy as np
import os,sys
from ase import Atoms, Atom
from ase.collections import methane
# from ase.db import connect
from ase.units import *
from reaction_tools import *
#
# Calculate entropy change along the reaction
#
argvs   = sys.argv
infile  = argvs[1]
outfile = "deltaS.txt"
f = open(outfile,"w")

# db = connect("methane_test.json")

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

rxn_num = get_number_of_reaction(infile)
#
# --- energy calculation ---
#
units   = create_units('2014')
amu     = units['_amu']
kbolt   = units['_k'] # kB in unit is not eV/K --> do not use
hplanck = units['_hplanck']
Nav     = units['_Nav']

type_for = ["gas"]*rxn_num
type_rev = ["gas"]*rxn_num

for irxn in range(rxn_num):
	reac_S = 0.0
	prod_S = 0.0
	#
	# reactants
	#
	for imol, mol in enumerate(r_ads[irxn]):
		mol  = mol[0]
		mol  = remove_side_and_flip(mol)

		if mol == 'surf' or mol == 'def':
			rxntype.append('surf')
		else:
			tmp = methane[mol]
			try:
				entropy = methane.data[mol]['molecular_entropy']
				print "molecular entropy for",mol," = ",entropy, "eV/K"
			except:
				entropy = 0.0

		reac_S += entropy

	#
	# products
	#
	for imol, mol in enumerate(p_ads[irxn]):
		mol  = mol[0]
		mol  = remove_side_and_flip(mol)

		if mol == 'surf' or mol == 'def':
			rxntype.append('surf')
		else:
			tmp = methane[mol]
			try:
				entropy = methane.data[mol]['molecular_entropy']
				print "molecular entropy for",mol," = ",entropy, "eV/K"
			except:
				entropy = 0.0

		prod_S += entropy

	deltaS = np.sum(prod_S) - np.sum(reac_S)

	# write to file
	f.write("{0:>16.8e}\n".format(deltaS))

