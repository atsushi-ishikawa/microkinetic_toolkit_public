import numpy as np
import os,sys
from reaction_tools import *

#
# calculate pre-exponential factor
#
reactionfile = "reaction.txt"
(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

rxn_num = get_number_of_reaction(reactionfile)

print r_ads, r_site, r_coef
print p_ads, p_site, p_coef

#
# --- energy calculation ---
#
from ase import Atoms, Atom
from ase.collections import methane

reac_A = np.array(rxn_num*[range(len(r_ads[0]))],dtype="f")
prod_A = np.array(rxn_num*[range(len(p_ads[0]))],dtype="f")

for irxn in range(rxn_num):
	print "---------", irxn, "------------"

	#
	# reactants
	#
	for imol, mol in enumerate(r_ads[irxn]):
		tmp  = methane[mol]
		mass = sum(tmp.get_masses())
		reac_A[irxn,imol] = mass
	#
	# products
	#
	for imol, mol in enumerate(p_ads[irxn]):
		tmp  = methane[mol]
		mass = sum(tmp.get_masses())
		prod_A[irxn,imol] = mass

