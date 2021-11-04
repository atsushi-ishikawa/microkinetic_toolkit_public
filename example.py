from reactions import Reaction, Reactions
from ase.build import fcc111

T = 800  # temperature [K]
P = 100  # total pressure [bar]

# read reactions from file and set Reactions
reactions = Reactions.from_csv("reactions.csv")

# define surface
surf = fcc111("Ni", size=[3, 3, 4], vacuum=10.0)

# adsorbate to surface?
reactions.adsorbate_on_surface(surf)

# get reaction energies for all the elementary reactions
deltaE = reactions.get_reaction_energies(method=emt)

# calculate rate constant from reaction energies
rateconst = reactions.get_rate_constants(deltaE=deltaE, T=T, P=P)

# solve rate equations
rate = reactions.solve_rate_equation(rateconstant=rateconst, T=T, P=P)

# draw graph
reactions.draw_graph(rate)

