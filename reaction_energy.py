import numpy as np
import sys

argvs = sys.argv
reac  = sys.argv[1]

print("calculating reaction %s") % reac

reac, prod = reac.split("-->")
reac = reac.replace(" ","")
prod = prod.replace(" ","")

reac = reac.split("+")
prod = prod.split("+")
reac_num = len(reac)
prod_num = len(prod)

print("reactant: %s, product: %s") % (reac,prod)
print reac_num, prod_num

r_ads = range(reac_num) ; r_site = range(reac_num)
for i,j in enumerate(reac):
	r_ads[i]  = j.split("_")[0]
	r_site[i] = j.split("_")[1]

p_ads = range(prod_num) ; p_site = range(prod_num)
for i,j in enumerate(prod):
 	p_ads[i]  = j.split("_")[0]
  	p_site[i] = j.split("_")[1]

print r_ads, r_site
print p_ads, p_site


