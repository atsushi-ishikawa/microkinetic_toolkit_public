def read_reactionfile(file):
	import os
	import matplotlib.pyplot as plt
	import networkx as nx
	import re

	os.system('grep -v "^#"    %s  > reaction2.txt' % file)
	os.system('grep -v "^\s*$" reaction2.txt > reaction3.txt')

	numlines = sum(1 for line in open("reaction3.txt"))

	f = open("reaction3.txt","r")

	lines = f.readlines()

	reac = range(numlines)
	rxn  = range(numlines)
	prod = range(numlines)

	for i,line in enumerate(lines):
		#text = line.replace("\n","").replace(">","").replace(" ","").split("--")
		text = line.replace("\n","").replace(">","").split("--")
		reac_tmp  = text[0]
		rxn_tmp   = text[1]
		prod_tmp  = text[2]

		reac[i] = re.split(" \+ ",reac_tmp) # for cations
		prod[i] = re.split(" \+ ",prod_tmp) # for cations

		reac = remove_space(reac)
		prod = remove_space(prod)

		rxn[i]  = reac[i][0] + "_" + rxn_tmp

	os.system('rm reaction2.txt reaction3.txt')

	return reac,rxn,prod

def remove_space(nested_list):
		li = nested_list
		for i,j in enumerate(li):
			tmp = j
			for ii,jj in enumerate(tmp):
				tmp[ii] = tmp[ii].strip()

			li[i] = tmp
		
		return li

