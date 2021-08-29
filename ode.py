import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pickle
from ode_equation import func

R = 8.314*1.0e-3  # kJ/mol/K
T = 500    # K

# get reaction energy
deltaE = []
f = open("deltaE.txt", "r")
for line in f:
	line = line.rstrip("\n")
	line = line.split()
	deltaE.append(float(line[1]))

Nrxn = len(deltaE)

# input by yourself
deltaS = [0, 0]

deltaE = np.array(deltaE)
deltaS = np.array(deltaS)

deltaG = deltaE - T*deltaS
eVtokJ = 96.487
deltaG = deltaG*eVtokJ

# get Ea
alpha = 1.0
beta  = 8.0

Ea = alpha*deltaE + beta

Afor = np.array([1.0e6]*Nrxn)
Kci  = np.exp(-deltaG/R/T)

sden  = 1.0e-14      # site density [mol/m^2]
w_cat = 100*1.0e-3   # catalyst weight [g]
area  = 1.0e4*w_cat*1.0e4  # surface area of catalyst [cm^2] = [cm^2/g] * [g]. Converted [cm^2]-->[m^2]

phi   = 0.5     # porosity
rho_b = 1.0e6   # density of catalyst [g/m^3]. typical is 1.0 g/cm^3 = 1.0*10^6 g/m^3
Vr    = w_cat/(rho_b*(1-phi))  # reactor volume [m^3], calculated from w_cat.

print("Ea:", Ea)
print("deltaG:", deltaG)
print("Kci:", Kci)

speciesfile = "species.pickle"
species = pickle.load(open(speciesfile, "rb"))

Ncomp = len(species)
Ngas  = len(list(filter(lambda x: "surf" not in x, species)))

t0, tf = 0, 1.0e-4
dt = 1.0e-8
t_span = (t0, tf)
t_eval = np.arange(t0, tf, dt)

print(species)
C0 = [0, 1, 1, 0, 0, 1]

soln = solve_ivp(fun=lambda t, C: func(t, C, Afor, Ea, Kci, T, sden, area, Vr, Ngas, Ncomp),
				 t_span=t_span, t_eval=t_eval, y0=C0,
				 rtol=1e-6, atol=1e-9, method="BDF")
print(soln.nfev, "evaluations requred.")

fig, (fig1, fig2) = plt.subplots(ncols=2, figsize=(10, 4))

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
