import numpy as np
import os,sys
#
# form reactant and product information
#
argvs = sys.argv
reac  = sys.argv[1]

print("calculating reaction %s") % reac

reac, prod = reac.split("-->")
reac = reac.replace(" ","")
prod = prod.replace(" ","")

reac = reac.split("+")
prod = prod.split("+")
reac_num = len(reac)
prod_num = len(prod)

#print("reactant: %s, product: %s") % (reac,prod)
#print reac_num, prod_num

r_ads = range(reac_num) ; r_site = range(reac_num)
for i,j in enumerate(reac):
	r_ads[i]  = j.split("_")[0]
	r_site[i] = j.split("_")[1]

p_ads = range(prod_num) ; p_site = range(prod_num)
for i,j in enumerate(prod):
 	p_ads[i]  =  j.split("_")[0]
  	p_site[i] = j.split("_")[1]

#print r_ads, r_site
#print p_ads, p_site

#
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.collections import methane
from ase.optimize import BFGS
from ase.vibrations import Vibrations
#
# now calculate reaction energy
# molecule's data should be stored in "methane.json"
#
method = "b3lyp"
basis  = "aug-cc-pvdz"
ZPE = True

reac_en = np.array(range(len(r_ads)),dtype="f")
prod_en = np.array(range(len(p_ads)),dtype="f")

for i,mol in enumerate(r_ads):
	tmp = methane[mol]
	natom = len(tmp.get_atomic_numbers())
	tmp.calc = Gaussian(method=method, basis=basis)
	opt = BFGS(tmp)
	opt.run(fmax=0.05)
	en  = tmp.get_potential_energy()
	if ZPE == True and natom != 1:
		vib = Vibrations(tmp)
		vib.run()
		hnu = vib.get_energies()
		zpe = vib.get_zero_point_energy()
		reac_en[i] = en + zpe
		os.system("rm vib.*")
	else:
		reac_en[i] = en

for i,mol in enumerate(p_ads):
	tmp = methane[mol]
	natom = len(tmp.get_atomic_numbers())
	tmp.calc = Gaussian(method=method, basis=basis)
	opt = BFGS(tmp)
	opt.run(fmax=0.05)
	en  = tmp.get_potential_energy()
	if ZPE == True and natom != 1:
		vib = Vibrations(tmp)
		vib.run()
		hnu = vib.get_energies()
		zpe = vib.get_zero_point_energy()
		prod_en[i] = en + zpe
		os.system("rm vib.*")
	else:
		prod_en[i] = en

deltaE = np.sum(prod_en) - np.sum(reac_en)

print deltaE

