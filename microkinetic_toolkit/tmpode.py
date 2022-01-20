import numpy as np

def func(t, c, kfor, Kc, T, sden, area, Vr, ngas, ncomp):

	krev = kfor / Kc
	theta = c[0:ncomp]
	theta = theta * sden
	rate = np.zeros(18)

	rate[0] = + 1.0*kfor[13]*c[6] - 1.0*krev[13]*c[0]*theta[17]*(area/Vr)  # CH4
	rate[1] = - 1.0*kfor[1]*c[1]*theta[17]*(area/Vr) + 1.0*krev[1]*c[9]*(area/Vr)  # CO2
	rate[2] = - 1.0*kfor[0]*c[2]*theta[17]**2*(area/Vr) + 1.0*krev[0]*c[14]**2*(area/Vr)  # H2
	rate[3] = + 1.0*kfor[16]*c[12] - 1.0*krev[16]*c[3]*theta[17]*(area/Vr)  # H2O
	rate[4] = + 1.0*kfor[10]*c[8]*c[14] - 1.0*krev[10]*c[4]*theta[17]*(area/Vr) - 1.0*kfor[11]*c[4]*c[14]*(area/Vr) + 1.0*krev[11]*c[5]*theta[17]*(area/Vr)  # CH2_surf
	rate[5] = + 1.0*kfor[11]*c[4]*c[14] - 1.0*krev[11]*c[5]*theta[17]*(area/Vr) - 1.0*kfor[12]*c[5]*c[14]*(area/Vr) + 1.0*krev[12]*c[6]*theta[17]*(area/Vr)  # CH3_surf
	rate[6] = + 1.0*kfor[12]*c[5]*c[14] - 1.0*krev[12]*c[6]*theta[17]*(area/Vr) - 1.0*kfor[13]*c[6]*(area/Vr) + 1.0*krev[13]*c[0]*theta[17]*(area/Vr)  # CH4_surf
	rate[7] = + 1.0*kfor[4]*c[13]*theta[17]*(area/Vr) - 1.0*krev[4]*c[7]*c[16]*(area/Vr) - 1.0*kfor[5]*c[7]*theta[17]*(area/Vr) + 1.0*krev[5]*c[10]*c[14]*(area/Vr) + 1.0*kfor[7]*c[10]*c[14]*(area/Vr) - 1.0*krev[7]*c[7]*theta[17]*(area/Vr) - 1.0*kfor[8]*c[7]*theta[17]*(area/Vr) + 1.0*krev[8]*c[11]*c[15]*(area/Vr)  # CHO_surf
	rate[8] = + 1.0*kfor[9]*c[11]*c[14] - 1.0*krev[9]*c[8]*theta[17]*(area/Vr) - 1.0*kfor[10]*c[8]*c[14]*(area/Vr) + 1.0*krev[10]*c[4]*theta[17]*(area/Vr)  # CH_surf
	rate[9] = + 1.0*kfor[1]*c[1]*theta[17]*(area/Vr) - 1.0*krev[1]*c[9]*(area/Vr) - 1.0*kfor[2]*c[9]*theta[17]*(area/Vr) + 1.0*krev[2]*c[10]*c[16]*(area/Vr) - 1.0*kfor[3]*c[9]*c[14]*(area/Vr) + 1.0*krev[3]*c[13]*theta[17]*(area/Vr)  # CO2_surf
	rate[10] = + 1.0*kfor[2]*c[9]*theta[17]*(area/Vr) - 1.0*krev[2]*c[10]*c[16]*(area/Vr) + 1.0*kfor[5]*c[7]*theta[17]*(area/Vr) - 1.0*krev[5]*c[10]*c[14]*(area/Vr) - 1.0*kfor[6]*c[10]*theta[17]*(area/Vr) + 1.0*krev[6]*c[11]*c[16]*(area/Vr) - 1.0*kfor[7]*c[10]*c[14]*(area/Vr) + 1.0*krev[7]*c[7]*theta[17]*(area/Vr)  # CO_surf
	rate[11] = + 1.0*kfor[6]*c[10]*theta[17]*(area/Vr) - 1.0*krev[6]*c[11]*c[16]*(area/Vr) + 1.0*kfor[8]*c[7]*theta[17]*(area/Vr) - 1.0*krev[8]*c[11]*c[15]*(area/Vr) - 1.0*kfor[9]*c[11]*c[14]*(area/Vr) + 1.0*krev[9]*c[8]*theta[17]*(area/Vr)  # C_surf
	rate[12] = + 1.0*kfor[15]*c[15]*c[14] - 1.0*krev[15]*c[12]*theta[17]*(area/Vr) - 1.0*kfor[16]*c[12]*(area/Vr) + 1.0*krev[16]*c[3]*theta[17]*(area/Vr)  # H2O_surf
	rate[13] = + 1.0*kfor[3]*c[9]*c[14] - 1.0*krev[3]*c[13]*theta[17]*(area/Vr) - 1.0*kfor[4]*c[13]*theta[17]*(area/Vr) + 1.0*krev[4]*c[7]*c[16]*(area/Vr)  # HCOO_surf
	rate[14] = + 2.0*kfor[0]*c[2]*theta[17]**2*(area/Vr) - 2.0*krev[0]*c[14]**2*(area/Vr) - 1.0*kfor[3]*c[9]*c[14]*(area/Vr) + 1.0*krev[3]*c[13]*theta[17]*(area/Vr) + 1.0*kfor[5]*c[7]*theta[17]*(area/Vr) - 1.0*krev[5]*c[10]*c[14]*(area/Vr) - 1.0*kfor[7]*c[10]*c[14]*(area/Vr) + 1.0*krev[7]*c[7]*theta[17]*(area/Vr) - 1.0*kfor[9]*c[11]*c[14]*(area/Vr) + 1.0*krev[9]*c[8]*theta[17]*(area/Vr) - 1.0*kfor[10]*c[8]*c[14]*(area/Vr) + 1.0*krev[10]*c[4]*theta[17]*(area/Vr) - 1.0*kfor[11]*c[4]*c[14]*(area/Vr) + 1.0*krev[11]*c[5]*theta[17]*(area/Vr) - 1.0*kfor[12]*c[5]*c[14]*(area/Vr) + 1.0*krev[12]*c[6]*theta[17]*(area/Vr) - 1.0*kfor[14]*c[16]*c[14]*(area/Vr) + 1.0*krev[14]*c[15]*theta[17]*(area/Vr) - 1.0*kfor[15]*c[15]*c[14]*(area/Vr) + 1.0*krev[15]*c[12]*theta[17]*(area/Vr)  # H_surf
	rate[15] = + 1.0*kfor[8]*c[7]*theta[17]*(area/Vr) - 1.0*krev[8]*c[11]*c[15]*(area/Vr) + 1.0*kfor[14]*c[16]*c[14]*(area/Vr) - 1.0*krev[14]*c[15]*theta[17]*(area/Vr) - 1.0*kfor[15]*c[15]*c[14]*(area/Vr) + 1.0*krev[15]*c[12]*theta[17]*(area/Vr)  # OH_surf
	rate[16] = + 1.0*kfor[2]*c[9]*theta[17]*(area/Vr) - 1.0*krev[2]*c[10]*c[16]*(area/Vr) + 1.0*kfor[4]*c[13]*theta[17]*(area/Vr) - 1.0*krev[4]*c[7]*c[16]*(area/Vr) + 1.0*kfor[6]*c[10]*theta[17]*(area/Vr) - 1.0*krev[6]*c[11]*c[16]*(area/Vr) - 1.0*kfor[14]*c[16]*c[14]*(area/Vr) + 1.0*krev[14]*c[15]*theta[17]*(area/Vr)  # O_surf
	rate[17] = -rate[4] -rate[5] -rate[6] -rate[7] -rate[8] -rate[9] -rate[10] -rate[11] -rate[12] -rate[13] -rate[14] -rate[15] -rate[16]  # surf

	# species --- 0 = CH4 1 = CO2 2 = H2 3 = H2O 4 = CH2_surf 5 = CH3_surf 6 = CH4_surf 7 = CHO_surf 8 = CH_surf 9 = CO2_surf 10 = CO_surf 11 = C_surf 12 = H2O_surf 13 = HCOO_surf 14 = H_surf 15 = OH_surf 16 = O_surf 17 = surf 
	if ncomp > ngas:
		rate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface
	return rate

