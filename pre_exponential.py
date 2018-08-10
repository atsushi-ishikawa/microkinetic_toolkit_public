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
T = 300.0 # temporary temperature

# Note on revserse reaction
# Pre-exponential factor for reverse reaction is generally not needed
# since reverse pre-exponential is calculated from forward one and equilibrium constant

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

rxn_num = get_number_of_reaction(infile)
#
# --- energy calculation ---
#
units   = create_units('2014')
amu     = units['_amu']
kbolt   = units['_k']
hplanck = units['_hplanck']
Nav     = units['_Nav']

type_for = ["gas"]*rxn_num

for irxn in range(rxn_num):
	#
	# reactants
	#
	mass_sum = 0; mass_prod = 1;
	rxntype = []

	sigmaAB = 1.0


	for imol, mol in enumerate(r_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)
		nmol = len(r_ads[irxn])

		if mol == 'surf' or mol == 'def':
			rxntype.append('surf')
		else:
			tmp = methane[mol]
			vol = methane.data[mol]['molecular_volume']
			sigma  = 3.0/4.0/np.pi * np.cbrt(vol) # molecular diagmeter (in Angstrom) calculated from radius
			sigma *= 10e-10 					  # Angstrom --> m
			coef   = r_coef[irxn][imol]
			sigma  = sigma**coef

			if nmol==1 and coef==1:
				sigmaAB  = sigma**2 * 10**0
			else:
				sigmaAB *= sigma

			mass = sum(tmp.get_masses())

			site = r_site[irxn][imol][0]
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

	# end loop over molecule

	if all(rxn == 'gas' for rxn in rxntype):
		#
		# gas reaction --- collision theory expression
		#
		type_for[irxn] = "gas"
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		fac_for  = sigmaAB * np.sqrt( kbolt*T/red_mass ) * Nav
		fac_for  = fac_for * 10**6 # [m^3] --> [cm^3]
		#
		# sqrt(kbolt*T/mass) [kg*m^2*s^-2*K^-1 * K * kg^-1]^1/2 = [m*s^-1]
		# A = [m^2] * [m*s^-1] * Nav = [m^3*s^-1]*[mol^-1] = [mol^-1*m^3*s^-1]
		# finally A in [mol^-1*cm^3*s^-1]
		#
	elif all(rxn == 'surf' for rxn in rxntype):
		if nmol == 1:
			#
			# desorption --- transition state theory
			#
			type_for[irxn] = "des"
			fac_for = kbolt/hplanck
		else:
			#
			# LH --- transition state theory
			#
			type_for[irxn] = "lh"
			fac_for = kbolt*T/hplanck # [s^-1]
	else:
		#
		# adsorption --- Hertz-Knudsen or Chemkin
		#
		type_for[irxn] = "ads"
		stick = 0.1 # rough !!!
		red_mass = mass_prod / mass_sum
		red_mass = red_mass*amu
		# --- Hertz-Knudsen
		# denom    = np.sqrt( 2.0*np.pi*red_mass*kbolt*T ) # [kg* kg*m^2*s^-2*K^-1 * K]^1/2 = [kg^2*m^2*s^-2]^1/2 = [kg*m*s^-1]
		# fac_for  = stick / denom # [kg^-1*m^-1*s]
		# fac_for *= 10**-5        # [g^-1*cm^-1*s]
		# --- chemkin
		sden  = 1.0*10**15 # [site/cm^2]
		sden /= Nav        # [site/cm^2]*[mol/site] = [mol/cm^2]
		sden *= 10**4      # [mol/cm^2] --> [mol/m^2]
		fac_for  = ( stick/sden ) * np.sqrt( kbolt*Nav*T / (2.0*np.pi*red_mass*Nav)) # [mol^-1*m^3*s^-1]
		fac_for *= 10**6 # [mol^-1*cm^3*s^-1]
		
	f.write("{0:>16.8e}\t{1:>6s}\n".format(fac_for, type_for[irxn]))

f.close()

