from ase.build import fcc111, surface, sort
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
import os

vacuum = 10.0

os.system('rm surf.db')

#surf =  fcc111(symbol="Cu", size=[3,3,4], a=3.6, vacuum=vacuum)
bulk = read("mgo.cif")
surf = surface(lattice=bulk, indices=(1,0,0), layers=1, vacuum=vacuum)
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

db = connect("surf.db")
db.write(surf, data=pos)

