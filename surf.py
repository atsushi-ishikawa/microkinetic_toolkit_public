from ase.build import fcc111, surface, sort, niggli_reduce
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
from ase import Atoms, Atom
import os
import numpy as np
from reaction_tools import *

vacuum = 10.0
doping = False

os.system('rm surf.db')

nlayer = 2
nrelax = 1

#cif_file = "mgo.cif"
cif_file = "La2O3.cif"
bulk = read(cif_file)
#niggli_reduce(bulk)
surf = surface(lattice=bulk, indices=(0,0,1), layers=nlayer, vacuum=vacuum) # step: (310) is good. nlayer=7, [1,2,1] might be good.

# 
if cif_file == "La2O3.cif":
	surf.rotate(180,'y', rotate_cell=False) # La2O3
	surf.wrap()
	surf.center(axis=2) # La2O3, only z-axis

surf = surf*[3,3,1]
#surf = sort(surf)
surf = sort_atoms_by_z(surf)

formula = surf.get_chemical_formula()

#
# doping e.g.) Mg by Li
#
if doping:
	symbols =  np.array(surf.get_chemical_symbols())
	#rep_atom = np.max( np.where(symbols == 'Mg') )
	rep_atom = 48  # 16 for layer=1, 48 for layer=2
	symbols[rep_atom] = 'Li'
	surf.set_chemical_symbols(symbols)

surf.translate([0,0,-vacuum+1])

#lattice = "fcc"
lattice = "hcp"
facet   = "001"

pos = {	'lattice' : lattice, 
		'facet'   : facet  ,
		'formula' : formula
      }

natoms = len(surf.get_atomic_numbers())
per_layer = natoms / nlayer / 2 # divide by 2 for oxides -- check
tag = np.ones(natoms, int)
if nrelax == 0:
	for i in range(natoms):
		tag[i] = 0
else:
	for i in range(natoms-1, natoms-nrelax*per_layer-1, -1):
		tag[i] = 0

surf.set_tags(tag)

db = connect("surf.db")
db.write(surf, data=pos)

