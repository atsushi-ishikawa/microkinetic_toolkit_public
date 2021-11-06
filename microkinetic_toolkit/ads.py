import os,sys
import re
from reaction_tools import *

argvs = sys.argv
infile = argvs[1]

f = open(infile,"r")

lines = f.readlines()

mol  = range(len(lines))
site = range(len(lines))

for i,line in enumerate(lines):
	aaa,bbb = line.replace("\n","").split(":")
	mol[i]  = remove_space(aaa)
	bbb     = remove_space(bbb)
	site[i] = bbb.split(",")

print mol
print site

'''
	reac = range(numlines)
	rxn  = range(numlines)
	prod = range(numlines)

	for i,line in enumerate(lines):
		text = line.replace("\n","").replace(">","").split("--")
		reac_tmp  = text[0]
		rxn_tmp   = text[1]
		prod_tmp  = text[2]

		reac[i] = re.split(" \+ ", reac_tmp) # for cations
		prod[i] = re.split(" \+ ", prod_tmp) # for cations

		reac[i] = remove_space(reac[i])
		prod[i] = remove_space(prod[i])

		rxn[i]  = reac[i][0] + "_" + rxn_tmp

	os.system('rm reaction2.txt reaction3.txt')

	return reac,rxn,prod

'''
