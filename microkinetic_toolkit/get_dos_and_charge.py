import numpy as np
import os
import sys
import re
import json
import math
from reaction_tools import *
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.collections import methane
from ase.optimize import BFGS, MDMin, FIRE
from ase.vibrations import Vibrations, Infrared
from ase.db import connect
from ase.io import read, write
from ase.build import add_adsorbate, sort
from ase.visualize import view
from ase.neb import NEB, NEBTools
argvs = sys.argv

calculator = "vasp"
calculator = calculator.lower()
#
# temprary database to avoid overlapping calculations
#
surface = True

if surface:
    dbfile = 'surf.db'
    dbfile = os.path.join(os.getcwd(), dbfile)
    db = connect(dbfile)
    surf = db.get_atoms(id=1)

# fix atoms
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 2])
surf.set_constraint(c)

# --- VASP ---
if "vasp" in calculator:
    xc = "rpbe"
    prec = "normal"
    encut = 400.0  # 213.0 or 400.0 or 500.0
    potim = 0.10
    nsw = 100
    nelmin = 5
    nelm = 40  # default:40
    ediff = 1.0e-5
    ediffg = -0.1
    kpts = [3, 3, 1]
    ismear = 0
    sigma = 0.10
    vacuum = 10.0  # for gas-phase molecules. surface vacuum is set by surf.py
    setups = None
    ivdw = 12
    ialgo = 48  # normal=38, veryfast=48
    npar = 12
    nsim = npar
    lwave = False
    #setups = {"O" : "_h"}

    method = xc
    basis = ""

    # DFT+U
    DFTU = False
    if DFTU:
        ldau = "true"
        ldautype = 2
        ldau_luj = {'La': {'L': 3, 'U': 3.5, 'J': 0.0},
                    'O': {'L': -1, 'U': 0.0, 'J': 0.0}}
        ialgo = 38
    # charge
    neutral = True

label = surf.get_chemical_formula()
magmom = surf.get_initial_magnetic_moments()
# spin switch
if int(math.ceil(sum(magmom))) != 0:
    print("has unpaired electron")
    ispin = 2
else:
    ispin = 1
#
# surface
#
if neutral:
    if DFTU:
        surf.calc = Vasp(label=label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
                         encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=False,
                         ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
                         ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj)  # DFT+U
    else:
        surf.calc = Vasp(label=label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
                         encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=False,
                         ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts)  # normal
else:
    nelect = get_number_of_valence_electrons(surf)
    nelect = nelect - charge

    if DFTU:
        surf.calc = Vasp(label=label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
                         encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=False,
                         ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
                         nelect=nelect, lmono="true", ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj)  # charge + DFT+U
    else:
        surf.calc = Vasp(label=label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
                         encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=False,
                         ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
                         nelect=nelect, lmono="true")  # charge

# optimize
surf.get_potential_energy()

# single point
prec_sp = "accurate"
ismear_sp = -5
kpts_sp = [5, 5, 1]

surf.calc = Vasp(label=label, prec=prec_sp, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
                 encut=encut, ismear=ismear_sp, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave,
                 ibrion=-1, potim=potim, nsw=0, ediff=ediff, ediffg=ediffg, kpts=kpts_sp,
                 lcharg=True, laechg=True)  # normal

surf.get_potential_energy()

os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHG WAVECAR REPORT >& /dev/null')

home = os.environ['HOME']
chgsum = os.path.join(home, "vasp/vtstscripts/vtstscripts-935/chgsum.pl")
bader = os.path.join(home, "vasp/bader/bader")

os.system('%s %s/AECCAR0 %s/AECCAR2' % (chgsum, label, label))
os.system('%s %s/CHGCAR -ref CHGCAR_sum' % (bader,  label))
