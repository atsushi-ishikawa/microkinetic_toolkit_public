import numpy as np

def func(t, c, kfor, Kc, T, sden, area, Vr, ngas, ncomp):

	R = 8.314e-3  # kJ/mol/K
	krev = kfor / Kc
	theta = c[0:ncomp]
	theta = theta * sden
	rate = np.zeros(5)

	rate[0] = - 1.0*kfor[0]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1]  # H2
	rate[1] = - 1.0*kfor[0]*c[1]*c[1]*c[1]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1]*c[1]*c[1]*c[1]  # H
	rate[2] = - 1.0*kfor[1]*c[3]*c[3] - 1.0*kfor[1]*c[3]*c[3]*c[3]  # O2
	rate[3] = - 1.0*kfor[1]*c[3]*c[3]*c[3]*c[3]*c[3] - 1.0*kfor[1]*c[3]*c[3]*c[3]*c[3]*c[3]*c[3]  # CO2
	rate[4] = - 1.0*kfor[0]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1] - 1.0*kfor[1]*c[3]*c[3] - 1.0*kfor[1]*c[3]*c[3]*c[3]  # CO

	# species --- 0 = H2 1 = H 2 = O2 3 = CO2 4 = CO 
	if ncomp > ngas:
		rate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface
	return rate

