import numpy as np
import os, sys, pickle, argparse
from ase import Atoms, Atom
from ase.collections import methane
from reaction_tools import *
#
# Calculate entropy change along the reaction
#
parser = argparse.ArgumentParser()
parser.add_argument("--reactionfile", required=True, help="file with elementary reactions")
argvs = parser.parse_args()
infile  = argvs.reactionfile

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

rxn_num = get_number_of_reaction(infile)

deltaS = np.zeros(rxn_num)
for irxn in range(rxn_num):
	reac_S = 0.0
	prod_S = 0.0
	#
	# reactants
	#
	for imol, mol in enumerate(r_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)
		site = r_site[irxn][imol][0]

		if site != 'gas' or mol == 'surf' or mol == 'def':
			# surface species
			entropy = 0.0
		else:
			tmp = methane[mol]
			try:
				entropy = methane.data[mol]['molecular_entropy']
			except:
				entropy = 0.0

		reac_S += entropy

	#
	# products
	#
	for imol, mol in enumerate(p_ads[irxn]):
		mol  = mol[0]
		mol  = remove_side_and_flip(mol)
		site = p_site[irxn][imol][0]

		if site != 'gas' or mol == 'surf' or mol == 'def':
			entropy = 0.0
		else:
			tmp = methane[mol]
			try:
				entropy = methane.data[mol]['molecular_entropy']
			except:
				entropy = 0.0

		prod_S += entropy

	deltaS[irxn] = np.sum(prod_S) - np.sum(reac_S)

pickle.dump(deltaS, open("deltaS.pickle", "wb"))
