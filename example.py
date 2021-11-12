from microkinetic_toolkit.reactions import Reactions
from ase.build import fcc111

T = 800  # temperature [K]
P = 100  # total pressure [bar]
ratio = {"H2": 1.1, "CO2": 1.0}  # partial pressure ratio
#
# read reactions from file and set Reactions
#
reactions = Reactions.from_csv("test.csv")

# define surface if you like
#surf = fcc111("Ni", size=[3, 3, 4], vacuum=10.0)

#
#  get reaction energies for all the elementary reactions
#
reactions.calculator = "emt"

# if you have pre-calculated ase.db
#reactions.ase_db = "ase.db"

deltaEs = reactions.get_reaction_energies()

# calculate rate constant from reaction energies
ks = reactions.get_rate_constants(deltaEs=deltaEs, T=T)

# do_microkinetics
reactions.do_microkinetics(deltaEs=deltaEs, ks=ks, T=T, P=P, ratio=ratio)

# draw graph
#reactions.draw_network(rate)
