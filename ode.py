import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def f(t, c, Afor, Ea, Kci):
	R = 8.314
	T = 400

	Ea = np.array(Ea)
	kfor = Afor * np.exp(-Ea / R / T)
	krev = kfor / Kci

	Rate = np.array([0]*5)
	sden = 1.0e2
	theta = c[-2:]*sden

	Rate[0] = -kfor[0]*c[0]*theta[1]
	Rate[1] = +kfor[1]*c[0]*theta[1] - kfor[2]*c[1]
	Rate[2] = +kfor[2]*c[1]

	Rate[3] = +kfor[0]*c[0]*theta[1]
	Rate[4] = -Rate[3]

	dcdt = Rate
	return dcdt

t0, tf = 0, 1.0e0
c0 = [1, 0, 0, 0, 1]
t_span = (t0, tf)
dt = 1.0e-3
t_eval = np.arange(t0, tf, dt)
Afor = np.array([1.0e6]*3)
Ea   = np.array([50e3, 50e3, 30e3])
Kci  = np.array([1.0, 1.0, 1.0])

soln = solve_ivp(fun=lambda t, c: f(t, c, Afor, Ea, Kci),
				 t_span=t_span, t_eval=t_eval, y0=c0,
				 method="LSODA", rtol=1e-3, atol=1e-6)
print(soln.nfev, "evaluations requred.")

fig, (fig1, fig2) = plt.subplots(ncols=2, figsize=(10, 4))
fig1.plot(soln.t, soln.y[0], label="[X]")
fig1.plot(soln.t, soln.y[1], label="[Y]")
fig1.plot(soln.t, soln.y[2], label="[Z]")
fig2.plot(soln.t, soln.y[3], label="theta_1")
fig2.plot(soln.t, soln.y[4], label="theta_vac")

fig1.set_xlabel("times /s")
fig1.set_ylabel("concentration /arb.units")
fig2.set_xlabel("times /s")
fig2.set_ylabel("concentration /arb.units")
fig1.legend()
fig2.legend()

plt.show()
