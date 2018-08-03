import numpy as np
import os,sys
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.db import connect
#
# calculate reaction energy
# molecule's data should be stored in "methane.json"
#
mol = sys.argv[1]
calculator = "gau" ; calculator = calculator.lower()
db = connect("methane_test.json")
#
add_volume  = True
add_entropy = True
#
###
## --- Gaussian ---
if "gau" in calculator:
	method = "b3lyp"
	basis  = "6-31G*"
else:
	raise RuntimeError('molecular volume only implemented in Gaussian')

#
# add volume and entropy information for all the molecules in database
#
Nmol = len(db)

for imol in range(Nmol+10): # add 10 for safety
	try:
		id = imol+1
		name = db.get(id=id)
		tmp  = db.get_atoms(id=id)
		print "adding volume and entropy",name
		magmom = tmp.get_initial_magnetic_moments()
		#
		# volume and molecular total entropy calculation
		#
		if add_volume:
			tmp.calc = Gaussian(method=method, basis=basis, volume='tight')
			tmp.get_potential_energy()
			vol = tmp.get_molecular_volume()

		if add_entropy:
			tmp.calc = Gaussian(method=method, basis=basis, force=None, freq='noraman')
			tmp.get_potential_energy()
			entropy = tmp.get_molecular_entropy()
		#
		# look for magmom
		#
		try:
			magmom = tmp.get_magnetic_moments()
		except:
			magmom = tmp.get_initial_magnetic_moments()

		tmp.set_initial_magnetic_moments(magmom)
		#
		# now write to database
		#
		if add_volume:
			if add_entropy:
				db.write(tmp, key_value_pairs={'name' : mol, 'molecular_volume' : vol, 'molecular_entropy' : entropy})
			else:
				db.write(tmp, key_value_pairs={'name' : mol, 'molecular_volume' : vol})

		oldid = id
		db.delete([oldid])

	except:
		print "has nothing --- go to next"


