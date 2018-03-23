from ase.build import fcc111, surface, sort
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
import os
import numpy as np

vacuum = 10.0

os.system('rm surf.db')

#surf =  fcc111(symbol="Cu", size=[3,3,4], a=3.6, vacuum=vacuum)
nlayer = 1
nrelax = 0

bulk = read("mgo.cif")
surf = surface(lattice=bulk, indices=(1,0,0), layers=nlayer, vacuum=vacuum)
surf = surf*[2,2,1]
surf = sort(surf)

formula = surf.get_chemical_formula()

surf.translate([0,0,-vacuum+1])

lattice = "fcc"
facet   = "100"

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

