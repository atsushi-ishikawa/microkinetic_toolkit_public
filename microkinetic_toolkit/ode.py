import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pickle
from ode_equation import func

np.set_printoptions(precision=3, linewidth=100)

# constants
R = 8.314*1.0e-3  # kJ/mol/K
eVtokJ = 96.487

# parameters
Pin = 1e5    # inlet pressure in Pascal
T = 700    # K
# volumetric flowrate [m^3/sec]. 1 [m^2/sec] = 1.0e6 [mL/sec] = 6.0e7 [mL/min]
v0 = 1e-5

sden = 1.0e-05  # site density [mol/m^2]
w_cat = 1.0e-3   # catalyst weight [kg]
area = 1000*w_cat  # surface area. [m^2/kg] (e.g. BET) * [kg] --> [m^2]

phi = 0.5     # porosity
# density of catalyst [kg/m^3]. typical is 1.0 g/cm^3 = 1.0*10^3 kg/m^3
rho_b = 1.0e3
Vr = (w_cat/rho_b)*(1-phi)  # reactor volume [m^3], calculated from w_cat.
# Vr = 0.01e-6  # [m^3]

alpha = 0.6
beta = 1.0

# read reaction energy
deltaE = []
f = open("deltaE.txt", "r")
for line in f:
    line = line.rstrip("\n")
    line = line.split()
    deltaE.append(float(line[1]))

Nrxn = len(deltaE)
deltaE = np.array(deltaE)
deltaE = deltaE*eVtokJ

# read species
speciesfile = "species.pickle"
species = pickle.load(open(speciesfile, "rb"))

Ncomp = len(species)
Ngas = len(list(filter(lambda x: "surf" not in x, species)))

# read entropy (in eV)
deltaS = pickle.load(open("deltaS.pickle", "rb"))
TdeltaS = T*deltaS
TdeltaS = TdeltaS*eVtokJ

# read pre-exponential (in [m, mol, s])
Afor, Afor_type = pickle.load(open("pre_exp.pickle", "rb"))
# multipy temperature and surface density
for i, itype in enumerate(Afor_type):
    if itype == "gas":
        Afor[i] *= np.sqrt(T)
    elif itype == "ads":
        Afor[i] *= np.sqrt(T)/sden
    elif itype == "lh" or itype == "des":
        Afor[i] *= (T/sden)

# calculate deltaG
deltaG = deltaE - TdeltaS

# get Ea
Ea = alpha*(deltaE/eVtokJ) + beta
Ea = Ea*eVtokJ

Kpi = np.exp(-deltaG/R/T)  # in pressure unit
# Kci = Kpi*(101325/R/T)     # convert to concentration unit
Kci = Kpi*(R*T/1)     # convert to concentration unit

# for i, itype in enumerate(Afor_type):
#	if itype == "ads" or itype == "lh" or itype == "des":
#		Kci[i] *= sden

tau = Vr/v0  # residence time [sec]
# tau = 1.0  # residence time [sec]

# output results here
print("Ea [kJ/mol]:", Ea)
print("deltaE [kJ/mol]:", deltaE)
print("TdeltaS [kJ/mol]:", TdeltaS)
print("deltaG [kJ/mol]:", deltaG)
print(
    "residence time [sec]: {0:5.3e}, GHSV [hour^-1]: {1:3d}".format(tau, int(60*60/tau)))

# now solve the ODE
t0, tf = 0, tau
dt = tf*1.0e-3
t_span = (t0, tf)
t_eval = np.arange(t0, tf, dt)

print("species:", species)
# C0 = PinPa / R*T
x0 = np.zeros(Ncomp)
x0[1] = 1.0  # CO2
x0[2] = 4.0  # H2

x0[-1] = 1.0  # vacancy
# density calculated from pressure. note that R is defined with kJ/mol/K.
C0 = Pin / (R*T*1e3)

# normalize x0 gas part
tot = np.sum(x0[:Ngas])
for i, j in enumerate(x0):
    if i <= Ngas:
        x0[i] = x0[i]/tot

C0 = x0*C0

soln = solve_ivp(fun=lambda t, C: func(t, C, Afor, Ea, Kci, T, sden, area, Vr, Ngas, Ncomp),
                 t_span=t_span, t_eval=t_eval, y0=C0,
                 rtol=1e-5, atol=1e-7, method="LSODA")  # method = {"BDF" | "Radau" | "LSODA"}
print(soln.nfev, "evaluations requred.")

fig, [fig1, fig2] = plt.subplots(ncols=2, figsize=(10, 4))

for i, isp in enumerate(species):
    if "surf" in isp:
        fig2.plot(soln.t, soln.y[i], label="theta{}".format(
            isp.replace("_", "").replace("surf", "")))
    else:
        fig1.plot(soln.t, soln.y[i], label="{}".format(isp))

fig1.set_xlabel("times /s")
fig1.set_ylabel("concentration /arb.units")
fig2.set_xlabel("times /s")
fig2.set_ylabel("concentration /arb.units")
fig1.legend()
fig2.legend()

plt.show()

# --- for graph plotting
make_rate_file = True
if make_rate_file:
    # variables
    pickle.dump((sden, area, Vr), open("variables.pickle", "wb"))

    # rate constant
    kfor = Afor * np.exp(-Ea/R/T)
    krev = kfor / Kci
    pickle.dump((kfor, krev), open("rateconst.pickle", "wb"))

#
# surface coverage: for time-independent, output the coverage at the last time step
#
tcov = []  # time-dependent coverage: tcov[species][time]
for i in range(Ncomp):
    tcov.append(soln.y[i])
tcov = np.array(tcov)

# time dependent coverage for graph
dT = 10
fac = 1.0e-6

coveragefile = "nodes.txt"
f = open(coveragefile, "wt")
f.write("      name      num     conc        time\n")  # header
# initial
for isp in range(Ngas):
    f.write("{0:16.14s}{1:03d}{2:12.4e}{3:12.4e}\n".format(
        species[isp], isp, tcov[isp][0], soln.t[0]))

for isp in range(Ngas, Ncomp):
    f.write("{0:16.14s}{1:03d}{2:12.4e}{3:12.4e}\n".format(
        species[isp], isp, tcov[isp][0]*fac, soln.t[0]))

# afterwards -- take averaged value
for it in range(dT, len(soln.y[0]), dT):
    for isp in range(Ngas):
        f.write("{0:16.14s}{1:03d}{2:12.4e}{3:12.4e}\n".format(
            species[isp], isp, tcov[isp][it], soln.t[it]))

    for isp in range(Ngas, Ncomp):
        f.write("{0:16.14s}{1:03d}{2:12.4e}{3:12.4e}\n".format(
            species[isp], isp, tcov[isp][it]*fac, soln.t[it]))

# add RXN nodes at last
for i in range(Nrxn):
    f.write("R{0:>03d}            {1:03d}{2:12.4e}{3:12.4e}\n".format(
        i, i, 1.0e-20, 0.0))

f.close()
