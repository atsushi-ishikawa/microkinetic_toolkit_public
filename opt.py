import numpy as np
import os, sys
from reaction_tools import *
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.calculators.nwchem import NWChem
from ase.calculators.emt import EMT
from ase.collections import methane
from ase.visualize import view
from ase.optimize import BFGS
from ase.db import connect
#
# calculate reaction energy
# molecule's data should be stored in "methane.json"
#
mol = sys.argv[1]
calculator = "gau"
calculator = calculator.lower()
db = connect("tmp.json")
###
## --- Gaussian ---
if "gau" in calculator:
	method = "b3lyp"
	basis  = "6-31G*"
## --- NWChem ---
elif "nw" in calculator:
	xc = "b3lyp"
	basis  = "6-31G*"
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
## --- EMT --- -> nothing to set

print("now optimize ", mol)
tmp    = methane[mol]
magmom = tmp.get_initial_magnetic_moments()
natom  = len(tmp.get_atomic_numbers())

if "gau" in calculator:
	tmp.calc = Gaussian(method=method, basis=basis)
	opt = BFGS(tmp)
	opt.run(fmax=0.05)
if "nw" in calculator:
	tmp.calc = NWChem(xc=xc, basis=basis, iterations=300)
	opt = BFGS(tmp)
	opt.run(fmax=0.05)
	#os.system('rm nwchem.*')
elif "vasp" in calculator:
	cell = [10.0, 10.0, 10.0]
	tmp.cell = cell
	tmp.calc = Vasp(prec=prec, xc=xc, ispin=2, encut=encut, ismear=0, istart=0,
					ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg,
					kpts=kpts)
	#tmp.center()
	tmp.get_potential_energy()
elif "emt" in calculator:
	tmp.calc = EMT()
	opt = BFGS(tmp)
	opt.run(fmax=0.05)

view(tmp)
db.write(tmp)
