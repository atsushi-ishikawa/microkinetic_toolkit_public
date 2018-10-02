from ase.build import fcc111, surface, sort, niggli_reduce
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
from ase import Atoms, Atom
from ase.constraints import FixAtoms
import os, sys
import numpy as np
import math
from reaction_tools import *
from ase.calculators.vasp import Vasp
from ase.units import J

argvs   = sys.argv
lattice = argvs[1]
facet   = argvs[2]

dbfile = "data.json"
db = connect(dbfile)

vacuum = 10.0
doping = False

nlayer = 4
nrelax = 1

cif_file = "mgo.cif"
#cif_file = "cao.cif"
#cif_file = "La2O3.cif"
#cif_file = "Ce2W3O12.cif"

name = cif_file.split(".")[0]

bulk = read(cif_file)

indices = []
for c in facet:
	indices.append(int(c))

surf = surface(lattice=bulk, indices=indices, layers=nlayer, vacuum=vacuum) # step: (310) is good. nlayer=7, [1,2,1] might be good.

if cif_file == "La2O3.cif":
	surf.rotate(180,'y', rotate_cell=False) # La2O3
	surf.wrap()
	surf.center(axis=2) # La2O3, only z-axis

surf = surf*[2,2,1]
#surf = surf*[1,2,1]
#surf = sort(surf)
surf = sort_atoms_by_z(surf)


formula = surf.get_chemical_formula()
#
# doping e.g.) Mg by Li
#
if doping:
	symbols  =  np.array(surf.get_chemical_symbols())
	rep_atom = 16 if nlayer == 1 else 48
	# MgO: 16 for layer=1, 48 for layer=2
	# CaO: 18 for layer=1
	symbols[rep_atom] = 'Li'
	surf.set_chemical_symbols(symbols)

#surf.translate([0,0,-vacuum+1])
surf.center()
#
# relaxation issues
#
natoms = len(surf.get_atomic_numbers())
per_layer = natoms / nlayer / 2 # divide by 2 for oxides -- check
#
# set tags: 2 = fixed layers, 1 = relaxing layers, 0 = adsorbates
#
# constraint will be set in reaction_energy.py
#
tag = np.full(natoms, 2, dtype=int)
if nrelax == 0:
	for i in range(natoms):
		tag[i] = 1
else:
	for i in range(0, nrelax*per_layer):
		tag[i] = 1
	for i in range(natoms-1, natoms-nrelax*per_layer-1, -1):
		tag[i] = 1

surf.set_tags(tag)

c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 2])
surf.set_constraint(c)
#view(surf) ; quit()
#
# do calculation
#
xc     = "pbesol"
prec   = "low"
encut  = 400.0 # 213.0 or 400.0 or 500.0
potim  = 0.15
nsw    = 10
nelmin = 5
nelm   = 20 # default:40
ediff  = 1.0e-4
ediffg = -0.3
kpts_blk = [3, 3, 3]
kpts_slb = [3, 3, 1]
ismear = 1
sigma  = 0.20
setups = None
ivdw   = 12
ialgo  = 48 # normal=38, veryfast=48
npar   = 12
nsim   = npar
lwave  = False
lcharg = True

bulk = bulk*[2,2,2]

label = name + "_" + lattice + "_" + facet + "_bulk"
bulk.calc = Vasp(label=label, prec=prec, xc=xc, ispin=2, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
				 encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
				 ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts_blk, isif=3)

label = name + "_" + lattice + "_" + facet + "_surf"
surf.calc = Vasp(label=label, prec=prec, xc=xc, ispin=2, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
				 encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
				 ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts_slb, isif=4)

Ebulk = bulk.get_potential_energy()
Ebulk_perN = Ebulk / len(bulk)
Eslab = surf.get_potential_energy()

cell = surf.get_cell_lengths_and_angles()
area = cell[0]*cell[1] * math.sin(math.radians(cell[3]))

print "bulk energy:",Ebulk
print "slab energy:",Eslab
print "area:",area
Esurf  = (Eslab - len(surf)*Ebulk_perN) / (2*area) # in ev/Ang^2
Esurf2 = Esurf/J*10**20 # in J/m^2
print("surface energy :%8.4f (eV/Ang^2),%8.4f (J/m^2):" % (Esurf, Esurf2))

data = {'lattice' : lattice,
		'facet'   : "[" + facet + "]",
		'Esurf'   : Esurf2	}
db.write(surf, data)

