from microkinetic_toolkit.reactions import Reactions
from microkinetic_toolkit.visualization.graph import ReactionsVisualizer
from ase.build import fcc111
import os
import pickle

# parameters
T = 900  # temperature [K]
P = 100  # total pressure [bar]
ratio = {"H2": 1.1, "CO2": 1.0}  # partial pressure ratio

# whether save files or not
deltaEs_pickle = "deltaEs.pickle"

# read reactions from file and set Reactions
reactions = Reactions.from_csv("test.csv")

# define surface if you like
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
# set parameters for kinetics
# BEP relatinoship (Ea = alpha*deltaE + beta)
bep_param = {"alpha": 0.84, "beta": 1.6}
sden  = 1.0e-5    # site density [mol/m^2]
v0    = 1.0e-6    # volumetric flowrate [m^3/sec]; 1.0e-6 [m^3/sec] = 1.0 [mL/sec]
area  = 1.0e3     # surface area (e.g. BET) [m^2/kg]. Typical value is 1.0e3.
wcat  = 1.0e-3    # catalyst weight [kg]
phi   = 0.5       # porocity
rho_b = 1.0e3     # density [kg/m^3]

reactions.set_kinetic_parameters(bep_param=bep_param, sden=sden, v0=v0, wcat=wcat, area=area, phi=phi, rho_b=rho_b)

# calculate rate constant from reaction energies
kfor = reactions.get_rate_constants(deltaEs=deltaEs, T=T)

# do_microkinetics
reactions.do_microkinetics(deltaEs=deltaEs, kfor=kfor, T=T, P=P, ratio=ratio)
#sens = reactions.get_sensitivity(deltaEs=deltaEs, kfor=kfor, T=T, P=P, ratio=ratio)
#reactions.get_rate_for_graph(T=T)

# draw graph
#graph = ReactionsVisualizer(reactions)
#graph.draw_network()
