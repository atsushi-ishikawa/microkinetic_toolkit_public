import numpy as np
import os,sys
from reaction_tools import *
#
# form reactant and product information
#
(reac, rxn, prod) = read_reactionfile("reaction.txt")

rxn_num = len(rxn)

r_ads  = [ [] for i in range(rxn_num) ]
p_ads  = [ [] for i in range(rxn_num) ]
r_site = [ [] for i in range(rxn_num) ]
p_site = [ [] for i in range(rxn_num) ]
r_coef = [ [] for i in range(rxn_num) ]
p_coef = [ [] for i in range(rxn_num) ]

for irxn, jrnx in enumerate(rxn):
	ireac = reac[irxn];     iprod = prod[irxn]
	ireac_num = len(ireac); iprod_num = len(iprod)

	r_ads[irxn]  = range(ireac_num)
	r_site[irxn] = range(ireac_num)
	r_coef[irxn] = range(ireac_num)

	for i,j in enumerate(ireac):
		if "*" in j:
			r_coef[irxn][i] = int( j.split("*")[0] )
			rest = j.split("*")[1]
		else:
			r_coef[irxn][i] = 1
			rest = j
		if "_" in rest:
			r_site[irxn][i] = rest.split("_")[1]
			r_ads[irxn][i]  = rest.split("_")[0]
		else:
			r_site[irxn][i] = "gas"
			r_ads[irxn][i]  = rest

	p_ads[irxn]  = range(iprod_num)
	p_site[irxn] = range(iprod_num)
	p_coef[irxn] = range(iprod_num)

	for i,j in enumerate(iprod):
		if "*" in j:
			p_coef[irxn][i] = int( j.split("*")[0] )
			rest = j.split("*")[1]
		else:
			p_coef[irxn][i] = 1
			rest = j
		if "_" in rest:
			p_site[irxn][i] = rest.split("_")[1]
			p_ads[irxn][i]  = rest.split("_")[0]
		else:
			p_site[irxn][i] = "gas"
			p_ads[irxn][i]  = rest

#print r_ads, r_site, r_coef
#print p_ads, p_site, p_coef


#
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.collections import methane
from ase.optimize import BFGS
from ase.vibrations import Vibrations
#
# now calculate reaction energy
# molecule's data should be stored in "methane.json"
#
calculator = "Gaussian" ; calculator = calculator.lower()
ZPE = False

## --- Gaussian ---
if "gau" in calculator:
	method = "b3lyp"
	basis  = "6-311G*"
## --- VASP ---
elif "vasp" in calculator:
	xc = "pbe"
	prec = "normal"
	encut = 400.0
	potim = 0.10
	nsw = 100
	ediff = 1.0e-4
	ediffg = -0.03
	kpts = [1, 1, 1]

for irxn in range(rxn_num):
	print "---------", irxn, "------------"

	reac_en = np.array(range(len(r_ads[irxn])),dtype="f")
	prod_en = np.array(range(len(p_ads[irxn])),dtype="f")

	for imol, mol in enumerate(r_ads[irxn]):
		tmp   = methane[mol]
		magmom = tmp.get_initial_magnetic_moments()
		natom = len(tmp.get_atomic_numbers())
		coef  = r_coef[irxn][imol]

		if "gau" in calculator:
			tmp.calc = Gaussian(method=method, basis=basis)
			opt = BFGS(tmp)
			opt.run(fmax=0.05)
		elif "vasp" in calculator:
			cell = [10.0, 10.0, 10.0]
			tmp.cell = cell
			tmp.calc = Vasp(prec=prec,xc=xc,ispin=2,encut=encut, ismear=0, istart=0,
					ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg,
					kpts=kpts )

		en  = tmp.get_potential_energy()

		if ZPE == True and natom != 1:
			vib = Vibrations(tmp)
			vib.run()
			hnu = vib.get_energies()
			zpe = vib.get_zero_point_energy()
			reac_en[imol] = en + zpe
			os.system("rm vib.*")
		else:
			reac_en[imol] = en

		reac_en[imol] = coef * reac_en[imol]

	for imol, mol in enumerate(p_ads[irxn]):
		tmp   = methane[mol]
		natom = len(tmp.get_atomic_numbers())
		coef  = p_coef[irxn][imol]

		if "gau" in calculator:
			tmp.calc = Gaussian(method=method, basis=basis)
			opt = BFGS(tmp)
			opt.run(fmax=0.05)
		elif "vasp" in calculator:
			cell = [10.0, 10.0, 10.0]
			tmp.cell = cell
			tmp.calc = Vasp(prec=prec,xc=xc,ispin=2,encut=encut, ismear=0, istart=0,
					ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg,
					kpts=kpts )

		en  = tmp.get_potential_energy()

		if ZPE == True and natom != 1:
			vib = Vibrations(tmp)
			vib.run()
			hnu = vib.get_energies()
			zpe = vib.get_zero_point_energy()
			prod_en[imol] = en + zpe
			os.system("rm vib.*")
		else:
			prod_en[imol] = en

		prod_en[imol] = coef * prod_en[imol]

	deltaE = np.sum(prod_en) - np.sum(reac_en)

	print r_ads[irxn], "--->", p_ads[irxn], " DeltaE =", deltaE

