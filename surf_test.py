from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.io import write

surf = fcc111(symbol="Pt", size=[3,3,4], a=3.8, vacuum=10.0)

write("surf.db", surf)
