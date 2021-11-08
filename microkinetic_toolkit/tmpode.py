import numpy as np

def func(t, c, kfor, Kc, T, sden, area, Vr, ngas, ncomp):

	R = 8.314e-3  # kJ/mol/K
	krev = kfor / Kc
	theta = c[0:ncomp]
	theta = theta * sden
	rate = np.zeros(5)

	rate[0] = - 1.0*kfor[1]*c[0]*c[0]*c[0]*c[0]*c[0] - 1.0*kfor[1]*c[0]*c[0]*c[0]*c[0]*c[0]*c[0]  # CO2
	rate[1] = - 1.0*kfor[0]*c[1]*c[1]*c[1]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1]*c[1]*c[1]*c[1]  # H
	rate[2] = - 1.0*kfor[1]*c[0]*c[0] - 1.0*kfor[1]*c[0]*c[0]*c[0]  # O2
	rate[3] = - 1.0*kfor[0]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1] - 1.0*kfor[1]*c[0]*c[0] - 1.0*kfor[1]*c[0]*c[0]*c[0]  # CO
	rate[4] = - 1.0*kfor[0]*c[1]*c[1] - 1.0*kfor[0]*c[1]*c[1]*c[1]  # H2

	# species --- 0 = CO2 1 = H 2 = O2 3 = CO 4 = H2 
	if ncomp > ngas:
		rate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface
	return rate

