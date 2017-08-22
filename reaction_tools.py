def read_reactionfile(file):
	import os
	import matplotlib.pyplot as plt
	import networkx as nx

	os.system('grep -v "^#"    %s  > reaction2.txt' % file)
	os.system('grep -v "^\s*$" reaction2.txt > reaction3.txt')

	numlines = sum(1 for line in open("reaction3.txt"))

	f = open("reaction3.txt","r")

	lines = f.readlines()

	reac = range(numlines)
	rxn  = range(numlines)
	prod = range(numlines)

	for i,line in enumerate(lines):
		text = line.replace("\n","").replace(">","").replace(" ","").split("--")
		reac_tmp  = text[0]
		rxn_tmp   = text[1]
		prod_tmp  = text[2]

		reac[i] = reac_tmp.split("+")
		rxn[i]  = reac[i][0] + "_" + rxn_tmp
		prod[i] = prod_tmp.split("+")

	os.system('rm reaction2.txt reaction3.txt')

	return reac,rxn,prod

