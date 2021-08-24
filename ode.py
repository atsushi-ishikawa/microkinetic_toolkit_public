import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ode_equation import func

t0, tf = 0, 1.0e0
C0 = [0, 0, 1, 0, 1]
t_span = (t0, tf)
dt = 1.0e-3
t_eval = np.arange(t0, tf, dt)
Afor = np.array([1.0e6]*2)
Ea   = np.array([50e3, 50e3])
Kci  = np.array([1.0, 1.0])
sden = 1.0e2

Ngas = 5
Ncomp = 5

soln = solve_ivp(fun=lambda t, C: func(t, C, Afor, Ea, Kci, sden, Ngas, Ncomp),
				 t_span=t_span, t_eval=t_eval, y0=C0,
				 method="LSODA", rtol=1e-3, atol=1e-6)
print(soln.nfev, "evaluations requred.")

fig, (fig1, fig2) = plt.subplots(ncols=2, figsize=(10, 4))
for i in range(0, Ngas):
	fig1.plot(soln.t, soln.y[i], label="[{}]".format(i))

for i in range(Ngas, Ncomp):
	fig2.plot(soln.t, soln.y[i], label="theta{}".format(i))

fig1.set_xlabel("times /s")
fig1.set_ylabel("concentration /arb.units")
fig2.set_xlabel("times /s")
fig2.set_ylabel("concentration /arb.units")
fig1.legend()
fig2.legend()

plt.show()
