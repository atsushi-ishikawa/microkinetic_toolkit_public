from ase.build import fcc111, surface, sort, niggli_reduce
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
from ase import Atoms, Atom
import os
import numpy as np
from reaction_tools import *

vacuum = 12.0
doping = True

os.system('rm surf.db')

nlayer = 2
nrelax = 1

cif_file = "mgo.cif"
#cif_file = "cao.cif"
#cif_file = "La2O3.cif"
#cif_file = "Ce2W3O12.cif"

bulk = read(cif_file)
#niggli_reduce(bulk)
surf = surface(lattice=bulk, indices=(1,0,0), layers=nlayer, vacuum=vacuum) # step: (310) is good. nlayer=7, [1,2,1] might be good.

lattice = "fcc"
facet   = "100"
#lattice = "hcp"
#facet   = "001"
#lattice = "sp15"
#facet = "010"

# 
if cif_file == "La2O3.cif":
	surf.rotate(180,'y', rotate_cell=False) # La2O3
	surf.wrap()
	surf.center(axis=2) # La2O3, only z-axis

surf = surf*[2,2,1]
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

surf.translate([0,0,-vacuum+1])
#
# information for JSON file
#
pos = {	'lattice' : lattice, 
		'facet'   : facet  ,
		'formula' : formula
      }
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
	for i in range(natoms-1, natoms-nrelax*per_layer-1, -1):
		tag[i] = 1

surf.set_tags(tag)

db = connect("surf.db")
db.write(surf, data=pos)

