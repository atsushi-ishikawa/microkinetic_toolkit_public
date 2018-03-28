import numpy as np
import os,sys
from ase import Atoms, Atom
from ase.collections import methane
from ase.units import *
from reaction_tools import *
#
# Calculate pre-exponential factor.
# Note that temperature should be multiplied at MATLAB
#
argvs   = sys.argv
infile  = argvs[1]
outfile = "pre_exp.txt"
f = open(outfile,"w")

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

sigmaAB = 1.0e-10 # collision radius [m] -- approximate value

reac_A = np.array(rxn_num*[range(len(r_ads[0]))],dtype="f")
prod_A = np.array(rxn_num*[range(len(p_ads[0]))],dtype="f")

type_for = ["gas"]*rxn_num
type_rev = ["gas"]*rxn_num

for irxn in range(rxn_num):
	#
	# reactants
	#
	mass_sum = 0; mass_prod = 1;
	rxntype = []
	for imol, mol in enumerate(r_ads[irxn]):
		mol = mol[0]
		nmol = len(r_ads[irxn])
		if mol == 'surf':
			rxntype.append('surf')
		else:
			tmp  = methane[mol]
			#site = r_site[irxn][imol]
			site = r_site[irxn][imol][0]

			mass = sum(tmp.get_masses())

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			mass_sum  = mass_sum  + mass
			mass_prod = mass_prod * mass

			if site=='gas':
				rxntype.append(site)
			else:
				rxntype.append('surf')

	if all(rxn == 'gas' for rxn in rxntype):
		# gas reaction
		type_for[irxn] = "gas"
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		fac_for = np.pi * sigmaAB**2 * kbolt**(-3.0/2.0) * np.sqrt(8.0/np.pi/red_mass)
	elif all(rxn == 'surf' for rxn in rxntype):
		if nmol == 1:
			# desorption
			type_for[irxn] = "des"
			red_mass = mass_prod / mass_sum
			red_mass = red_mass*amu
			denom = np.sqrt( 2.0*np.pi*red_mass*kbolt )
			fac_for = 1.0 / denom
		else:
			# LH
			type_for[irxn] = "lh"
			fac_for = kbolt/hplanck
	else:
		# adsorption
		type_for[irxn] = "ads"
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		denom = np.sqrt( 2.0*np.pi*red_mass*kbolt )
		fac_for = 1.0 / denom

	#
	# products
	#
	mass_sum = 0; mass_prod = 1;
	rxntype = []
	for imol, mol in enumerate(p_ads[irxn]):
		mol = mol[0]
		nmol = len(p_ads[irxn])
		if mol == 'surf':
			rxntype.append('surf')
		else:
			tmp  = methane[mol]
			#site = p_site[irxn][imol]
			site = p_site[irxn][imol][0]

			mass = sum(tmp.get_masses())

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			mass_sum  = mass_sum  + mass
			mass_prod = mass_prod * mass

			if site=='gas':
				rxntype.append(site)
			else:
				rxntype.append('surf')

	if all(rxn == 'gas' for rxn in rxntype):
		# gas reaction
		type_rev[irxn] = "gas"
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		fac_rev = np.pi * sigmaAB**2 * kbolt**(-3.0/2.0) * np.sqrt(8.0/np.pi/red_mass)
	elif all(rxn == 'surf' for rxn in rxntype):
		if nmol == 1:
			# desorption
			type_rev[irxn] = "des"
			red_mass = mass_prod / mass_sum
			red_mass = red_mass*amu
			denom = np.sqrt( 2.0*np.pi*red_mass*kbolt )
			fac_rev = 1.0 / denom
		else:
			# LH
			type_rev[irxn] = "lh"
			fac_rev = kbolt/hplanck
	else:
		# adsorption
		type_rev[irxn] = "ads"
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		denom = np.sqrt( 2.0*np.pi*red_mass*kbolt )
		fac_rev = 1.0 / denom

	# write to file
	f.write("{0:>16.8e}\t{1:>16.8e}\t{2:6s}\t{3:6s}\n".format(fac_for,fac_rev,type_for[irxn],type_rev[irxn]))

