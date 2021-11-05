from reactions import Reaction, Reactions
from ase.build import fcc111

T = 800  # temperature [K]
P = 100  # total pressure [bar]

# read reactions from file and set Reactions
reactions = Reactions.from_csv("test.csv")

# define surface
surf = fcc111("Ni", size=[3, 3, 4], vacuum=10.0)

# get reaction energies for all the elementary reactions
deltaEs = reactions.get_reaction_energies(surface=surf, method="emt")

# calculate rate constant from reaction energies
rateconsts = reactions.get_rate_constants(deltaEs=deltaEs, T=T, P=P)

# do_microkinetics
reactions.do_microkinetics(rateconstants=rateconsts, T=T, P=P)

# draw graph
#reactions.draw_network(rate)
