from ase import Atoms,Atom
from ase.optimize.lbfgs import LBFGS
from ase.calculators.gaussian import Gaussian
from ase.visualize import view

h2 = Atoms("H2",[(0,0,0),(0,0,1.0)])
h2.calc = Gaussian(method="b3lyp",label="h2",basis="sto-3g")

opt = LBFGS(h2)
opt.run(fmax=0.05)

epot = h2.get_potential_energy()

view(h2)

# print epot

