from microkinetic_toolkit.reactions import Reactions
from ase.build import fcc111

# parameters
T = 900  # temperature [K]
P = 100  # total pressure [bar]
ratio = {"H2": 1.1, "CO2": 1.0}  # partial pressure ratio

## read reactions from file and set Reactions
reactions = Reactions.from_csv("test.csv")

## define surface if you like
reactions.surface = fcc111("Ni", size=[3, 3, 4], vacuum=10.0, periodic=True)

I_have_precalculated_database = False
if I_have_precalculated_database:
	reactions.ase_db = "ase.db"

# reaction energy evaluation
reactions.calculator = "emt"
deltaEs = reactions.get_reaction_energies()

# microkinetic analysis
## set parameters for kinetics
bep_param = {"alpha": 0.6, "beta": 1.0}  # BEP relatinoship (Ea = alpha*deltaE + beta)
sden = 1.0e-4  # site density
v0 = 1.0e-3  # volumetric flowrate [m^3/sec]
wcat = 100.0e-3  # catalyst weight [mg]
phi = 0.5  # porocity
rho_b = 1.0e3  # density

reactions.set_kinetic_parameters(bep_param=bep_param, sden=sden, v0=v0, wcat=wcat, phi=phi, rho_b=rho_b)

## calculate rate constant from reaction energies
ks = reactions.get_rate_constants(deltaEs=deltaEs, T=T)

## do_microkinetics
reactions.do_microkinetics(deltaEs=deltaEs, ks=ks, T=T, P=P, ratio=ratio)

## draw graph
#reactions.draw_network(rate)
