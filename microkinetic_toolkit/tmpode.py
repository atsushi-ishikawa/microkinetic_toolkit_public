import numpy as np

def func(t, c, kfor, Kc, T, sden, area, Vr, ngas, ncomp):

	krev = kfor / Kc
	theta = c[0:ncomp]
	theta = theta * sden
	rate = np.zeros(14)

	rate[0] = + 1.0*kfor[5]*theta[5]*(area/Vr) - 1.0*krev[5]*c[0]*theta[13]*(area/Vr)  # CH4
	rate[1] = - 1.0*kfor[1]*c[1]*theta[13]*(area/Vr) + 1.0*krev[1]*theta[6]*(area/Vr)  # CO2
	rate[2] = - 1.0*kfor[0]*c[2]*theta[13]**2*(area/Vr) + 1.0*krev[0]*theta[10]**2*(area/Vr)  # H2
	rate[3] = + 1.0*kfor[8]*theta[9]*(area/Vr) - 1.0*krev[8]*c[3]*theta[13]*(area/Vr)  # H2O
	rate[4] = - 1.0*kfor[4]*theta[4]*theta[10]*(area/Vr) + 1.0*krev[4]*theta[5]*theta[13]*(area/Vr)  # CH3_surf
	rate[5] = + 1.0*kfor[4]*theta[4]*theta[10]*(area/Vr) - 1.0*krev[4]*theta[5]*theta[13]*(area/Vr) - 1.0*kfor[5]*theta[5]*(area/Vr) + 1.0*krev[5]*c[0]*theta[13]*(area/Vr)  # CH4_surf
	rate[6] = + 1.0*kfor[1]*c[1]*theta[13]*(area/Vr) - 1.0*krev[1]*theta[6]*(area/Vr) - 1.0*kfor[2]*theta[6]*theta[13]*(area/Vr) + 1.0*krev[2]*theta[7]*theta[12]*(area/Vr)  # CO2_surf
	rate[7] = + 1.0*kfor[2]*theta[6]*theta[13]*(area/Vr) - 1.0*krev[2]*theta[7]*theta[12]*(area/Vr) - 1.0*kfor[3]*theta[7]*theta[13]*(area/Vr) + 1.0*krev[3]*theta[8]*theta[12]*(area/Vr)  # CO_surf
	rate[8] = + 1.0*kfor[3]*theta[7]*theta[13]*(area/Vr) - 1.0*krev[3]*theta[8]*theta[12]*(area/Vr)  # C_surf
	rate[9] = + 1.0*kfor[7]*theta[11]*theta[10]*(area/Vr) - 1.0*krev[7]*theta[9]*theta[13]*(area/Vr) - 1.0*kfor[8]*theta[9]*(area/Vr) + 1.0*krev[8]*c[3]*theta[13]*(area/Vr)  # H2O_surf
	rate[10] = + 2.0*kfor[0]*c[2]*theta[13]**2*(area/Vr) - 2.0*krev[0]*theta[10]**2*(area/Vr) - 1.0*kfor[4]*theta[4]*theta[10]*(area/Vr) + 1.0*krev[4]*theta[5]*theta[13]*(area/Vr) - 1.0*kfor[6]*theta[12]*theta[10]*(area/Vr) + 1.0*krev[6]*theta[11]*theta[13]*(area/Vr) - 1.0*kfor[7]*theta[11]*theta[10]*(area/Vr) + 1.0*krev[7]*theta[9]*theta[13]*(area/Vr)  # H_surf
	rate[11] = + 1.0*kfor[6]*theta[12]*theta[10]*(area/Vr) - 1.0*krev[6]*theta[11]*theta[13]*(area/Vr) - 1.0*kfor[7]*theta[11]*theta[10]*(area/Vr) + 1.0*krev[7]*theta[9]*theta[13]*(area/Vr)  # OH_surf
	rate[12] = + 1.0*kfor[2]*theta[6]*theta[13]*(area/Vr) - 1.0*krev[2]*theta[7]*theta[12]*(area/Vr) + 1.0*kfor[3]*theta[7]*theta[13]*(area/Vr) - 1.0*krev[3]*theta[8]*theta[12]*(area/Vr) - 1.0*kfor[6]*theta[12]*theta[10]*(area/Vr) + 1.0*krev[6]*theta[11]*theta[13]*(area/Vr)  # O_surf
	rate[13] = -rate[4] -rate[5] -rate[6] -rate[7] -rate[8] -rate[9] -rate[10] -rate[11] -rate[12]  # surf

	# species --- 0 = CH4 1 = CO2 2 = H2 3 = H2O 4 = CH3_surf 5 = CH4_surf 6 = CO2_surf 7 = CO_surf 8 = C_surf 9 = H2O_surf 10 = H_surf 11 = OH_surf 12 = O_surf 13 = surf 
	if ncomp > ngas:
		rate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface
	return rate

