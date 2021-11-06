import numpy as np

def func(t, c, k, Kc, T, sden, area, Vr, Ngas, Ncomp):
			R  = 8.314e-3  # kJ/mol/K
		
			krev = kfor / Kc
		
			theta = c[0:Ncomp]
			theta = theta * sden
		
