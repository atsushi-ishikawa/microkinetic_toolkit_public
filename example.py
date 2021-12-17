from microkinetic_toolkit.reactions import Reactions
from ase.build import fcc111
import os
import pickle

# parameters
T = 800  # temperature [K]
P = 100  # total pressure [bar]
ratio = {"H2": 1.1, "CO2": 1.0}  # partial pressure ratio

# whether save files or not
deltaEs_pickle = "deltaEs.pickle"

## read reactions from file and set Reactions
reactions = Reactions.from_csv("test.csv")

## define surface if you like
reactions.surface = fcc111("Ni", size=[3, 3, 4], vacuum=10.0, periodic=True)

I_have_precalculated_database = False
if I_have_precalculated_database:
	reactions.ase_db = "ase.db"

# reaction energy evaluation
reactions.calculator = "emt"
if os.path.exists(deltaEs_pickle):
	deltaEs = pickle.load(open(deltaEs_pickle, "rb"))
else:
	deltaEs = reactions.get_reaction_energies()
	with open(deltaEs_pickle, "wb") as f:
		pickle.dump(deltaEs, f)

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
kfor = reactions.get_rate_constants(deltaEs=deltaEs, T=T)

## do_microkinetics
reactions.do_microkinetics(deltaEs=deltaEs, kfor=kfor, T=T, P=P, ratio=ratio)

reactions.get_rate_for_graph(T=T)

## draw graph
reactions.draw_network()
