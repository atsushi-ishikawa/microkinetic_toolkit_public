from ase.build import fcc111, surface
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read
import os

vacuum = 10.0

os.system('rm surf.db')

#surf =  fcc111(symbol="Cu", size=[3,3,4], a=3.6, vacuum=vacuum)
bulk = read("mgo.cif")
surf = surface(lattice=bulk, indices=(1,0,0), layers=3, vacuum=vacuum)
surf = surf*[3,3,1]

surf.translate([0,0,-vacuum])

lattice = "fcc"
facet   = "100"

pos = {	'lattice' : lattice, 
		'facet'   : facet }

db = connect("surf.db")
db.write(surf, data=pos)

