import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pickle
from ode_equation import func

R   = 8.314*1.0e-3  # kJ/mol/K
Pin = 1e5  # inlet pressure in Pascal
T   = 900    # K

# get reaction energy
deltaE = []
f = open("deltaE.txt", "r")
for line in f:
	line = line.rstrip("\n")
	line = line.split()
	deltaE.append(float(line[1]))

Nrxn = len(deltaE)

# input by yourself (in eV)
deltaS = np.zeros(Nrxn)
deltaS[0] = -1.0e-3
deltaS[1] = -1.0e-3
deltaS[7] = 1.0e-3
deltaS[9] = 1.0e-3

deltaE = np.array(deltaE)
deltaS = np.array(deltaS)

# get detaG
deltaG = deltaE - T*deltaS
eVtokJ = 96.487
deltaG = deltaG*eVtokJ

# get Ea
alpha = 0.8
beta  = 1.2
Ea = alpha*deltaE + beta
Ea = Ea*eVtokJ

Afor = np.array([1.0e08]*Nrxn)  # in [m, mol, s]
Kci  = np.exp(-deltaG/R/T)

sden  = 1.0e-05  # site density [mol/m^2]
w_cat = 0.1e-3   # catalyst weight [kg]
area  = 1000*w_cat  # surface area. [m^2/kg] (e.g. BET) * [kg] --> [m^2]

phi   = 0.5     # porosity
rho_b = 1.0e3   # density of catalyst [kg/m^3]. typical is 1.0 g/cm^3 = 1.0*10^3 kg/m^3
Vr    = (w_cat/rho_b)*(1-phi)  # reactor volume [m^3], calculated from w_cat.
#Vr = 0.1e-6  # [m^3]

print("Ea:", Ea)
print("deltaG:", deltaG)
print("Kci:", Kci)

speciesfile = "species.pickle"
species = pickle.load(open(speciesfile, "rb"))

Ncomp = len(species)
Ngas  = len(list(filter(lambda x: "surf" not in x, species)))

t0, tf = 0, 1.0e1
dt = tf*1.0e-3
t_span = (t0, tf)
t_eval = np.arange(t0, tf, dt)

print(species)
# C0 = PinPa / R*T
x0 = np.zeros(Ncomp)
x0[1]  = 1.0  # CO2
x0[2]  = 1.0  # H2

x0[-1] = 1.0  # vacancy
C0 = Pin / (R*T*1e3)  # density calculated from pressure. note that R is defined with kJ/mol/K.
tot = np.sum(x0[:Ngas])

# normalize x0 gas part
for i, j in enumerate(x0):
	if i <= Ngas:
		x0[i] = x0[i]*C0/tot
	#else:
	#	x0[i] = x0[i]*sden

C0 = x0

soln = solve_ivp(fun=lambda t, C: func(t, C, Afor, Ea, Kci, T, sden, area, Vr, Ngas, Ncomp),
				 t_span=t_span, t_eval=t_eval, y0=C0,
				 rtol=1e-5, atol=1e-8, method="BDF")
print(soln.nfev, "evaluations requred.")

fig, [fig1, fig2] = plt.subplots(ncols=2, figsize=(10, 4))

for i, ispecies in enumerate(species):
	if "surf" in ispecies:
		fig2.plot(soln.t, soln.y[i], label="theta{}".format(ispecies.replace("_", "").replace("surf", "")))
	else:
		fig1.plot(soln.t, soln.y[i], label="{}".format(ispecies))

fig1.set_xlabel("times /s")
fig1.set_ylabel("concentration /arb.units")
fig2.set_xlabel("times /s")
fig2.set_ylabel("concentration /arb.units")
fig1.legend()
fig2.legend()

plt.show()
