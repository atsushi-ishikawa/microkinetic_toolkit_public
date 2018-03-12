import numpy as np
import os,sys
from reaction_tools import *
#
# calculate pre-exponential factor
#
#reactionfile = "reaction.txt"
argvs  = sys.argv
infile = argvs[1]
(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

rxn_num = get_number_of_reaction(infile)
#
# --- energy calculation ---
#
from ase import Atoms, Atom
from ase.collections import methane
from ase.units import kB

reac_A = np.array(rxn_num*[range(len(r_ads[0]))],dtype="f")
prod_A = np.array(rxn_num*[range(len(p_ads[0]))],dtype="f")

Temp = 300.0

for irxn in range(rxn_num):
	print "---------", irxn, "------------"
	#
	# reactants
	#
	mass_sum = 0; mass_prod = 1;
	for imol, mol in enumerate(r_ads[irxn]):
		print " === ", imol, "==="
		tmp  = methane[mol]
		site = r_site[irxn][imol]

		mass = sum(tmp.get_masses())
		if site=='gas':
			mass_sum  = mass_sum  + mass
			mass_prod = mass_prod * mass

		red_mass = mass_prod / mass_sum
		print red_mass
		print np.sqrt( 2.0*red_mass*np.pi * kB * Temp)
'''
	#
	# products
	#
	for imol, mol in enumerate(p_ads[irxn]):
		tmp  = methane[mol]
		mass = sum(tmp.get_masses())
		prod_A[irxn,imol] = mass

'''
