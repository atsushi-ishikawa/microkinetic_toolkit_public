import numpy as np
import os,sys
from reaction_tools import *
#
# calculate pre-exponential factor
#
#reactionfile = "reaction.txt"
argvs   = sys.argv
infile  = argvs[1]
outfile = "pre_exp.txt"
f = open(outfile,"w")

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

rxn_num = get_number_of_reaction(infile)
#
# --- energy calculation ---
#
from ase import Atoms, Atom
from ase.collections import methane
from ase.units import *

units = create_units('2014')
amu   = units['_amu']

reac_A = np.array(rxn_num*[range(len(r_ads[0]))],dtype="f")
prod_A = np.array(rxn_num*[range(len(p_ads[0]))],dtype="f")

rxn_type = np.array(rxn_num)
quit()

for irxn in range(rxn_num):
	#
	# reactants
	#
	mass_sum = 0; mass_prod = 1;
	for imol, mol in enumerate(r_ads[irxn]):
		tmp  = methane[mol]
		site = r_site[irxn][imol]

		mass = sum(tmp.get_masses())

		try:
			site,site_pos = site.split(".")
		except:
			site_pos = 'x1y1'

		if site=='gas':
			mass_sum  = mass_sum  + mass
			mass_prod = mass_prod * mass

		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		denom = np.sqrt( 2.0*np.pi*red_mass*kB )
		fac_for = 1.0 / denom
	#
	# products
	#
	mass_sum = 0; mass_prod = 1;
	for imol, mol in enumerate(p_ads[irxn]):
		tmp  = methane[mol]
		site = p_site[irxn][imol]

		mass = sum(tmp.get_masses())

		try:
			site,site_pos = site.split(".")
		except:
			site_pos = 'x1y1'

		if site=='gas':
			mass_sum  = mass_sum  + mass
			mass_prod = mass_prod * mass

		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		denom = np.sqrt( 2.0*np.pi*red_mass*kB )
		fac_rev = 1.0 / denom

	f.write("{0:>16.8e}\t{1:>16.8e}\n".format(fac_for,fac_rev))

