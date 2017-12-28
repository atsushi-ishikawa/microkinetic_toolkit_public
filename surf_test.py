from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.db import connect

surf    =  fcc111(symbol="Cu", size=[3,3,4], a=3.6, vacuum=10.0)

lattice = "fcc"
facet   = "111"

pos = {	'lattice' : lattice, 
		'facet'   : facet }

db = connect("surf.db")
db.write(surf, data=pos)

