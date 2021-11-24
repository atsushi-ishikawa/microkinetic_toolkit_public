from microkinetic_toolkit.reactions import Reactions
from ase.build import fcc111

# preparation
T = 800  # temperature [K]
P = 100  # total pressure [bar]
ratio = {"H2": 1.1, "CO2": 1.0}  # partial pressure ratio

## read reactions from file and set Reactions
reactions = Reactions.from_csv("test.csv")

## define surface if you like
#surf = fcc111("Ni", size=[3, 3, 4], vacuum=10.0)

## if you have pre-calculated ase.db
#reactions.ase_db = "ase.db"

# reaction energy evaluation
reactions.calculator = "emt"
deltaEs = reactions.get_reaction_energies()

# microkinetic analysis
## set parameters for kinetics
reactions.set_kinetic_parameters(alpha=0.5, beta=2.5, sden=1.0e-5, v0=1.0e-5, wcat=1.0e-3, phi=0.5, rho_b=1.0e3)

## calculate rate constant from reaction energies
ks = reactions.get_rate_constants(deltaEs=deltaEs, T=T)

## do_microkinetics
reactions.do_microkinetics(deltaEs=deltaEs, ks=ks, T=T, P=P, ratio=ratio)

## draw graph
#reactions.draw_network(rate)
