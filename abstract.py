import numpy as np
import os,sys
from reaction_tools import *
#
# form reactant and product information
#
(reac, rxn, prod) = read_reactionfile("reaction.txt")

rxn_num = len(rxn)

r_ads  = [ [] for i in range(rxn_num) ]
p_ads  = [ [] for i in range(rxn_num) ]
r_site = [ [] for i in range(rxn_num) ]
p_site = [ [] for i in range(rxn_num) ]
r_coef = [ [] for i in range(rxn_num) ]
p_coef = [ [] for i in range(rxn_num) ]

for irxn, jrnx in enumerate(rxn):
	ireac = reac[irxn];     iprod = prod[irxn]
	ireac_num = len(ireac); iprod_num = len(iprod)

	r_ads[irxn]  = range(ireac_num)
	r_site[irxn] = range(ireac_num)
	r_coef[irxn] = range(ireac_num)

	for i,j in enumerate(ireac):
		if "*" in j:
			r_coef[irxn][i] = int( j.split("*")[0] )
			rest = j.split("*")[1]
		else:
			r_coef[irxn][i] = 1
			rest = j
		if "_" in rest:
			r_site[irxn][i] = rest.split("_")[1]
			r_ads[irxn][i]  = rest.split("_")[0]
		else:
			r_site[irxn][i] = "gas"
			r_ads[irxn][i]  = rest

	p_ads[irxn]  = range(iprod_num)
	p_site[irxn] = range(iprod_num)
	p_coef[irxn] = range(iprod_num)

	for i,j in enumerate(iprod):
		if "*" in j:
			p_coef[irxn][i] = int( j.split("*")[0] )
			rest = j.split("*")[1]
		else:
			p_coef[irxn][i] = 1
			rest = j
		if "_" in rest:
			p_site[irxn][i] = rest.split("_")[1]
			p_ads[irxn][i]  = rest.split("_")[0]
		else:
			p_site[irxn][i] = "gas"
			p_ads[irxn][i]  = rest

print r_ads, r_site, r_coef
print p_ads, p_site, p_coef

