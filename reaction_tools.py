def read_reactionfile(file):
	import os
	import networkx as nx
	import re

	os.system('grep -v "^#"    %s  > reaction2.txt' % file)
	os.system('grep -v "^\s*$" reaction2.txt > reaction3.txt')

	numlines = sum(1 for line in open("reaction3.txt"))

	f = open(file,"r")

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

		reac[i] = re.split(" \+ ", reac_tmp) # for cations
		prod[i] = re.split(" \+ ", prod_tmp) # for cations

		reac[i] = remove_space(reac[i])
		prod[i] = remove_space(prod[i])

		rxn[i]  = reac[i][0] + "_" + rxn_tmp

	os.system('rm reaction2.txt reaction3.txt')

	return reac,rxn,prod

def read_speciesfile(file):
	import os
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

def remove_space(obj):
		newobj = [0]*len(obj)
		if isinstance(obj, str):
			#
			# string
			#
			newobj = obj.replace(" ","")
		elif isinstance(obj, list):
			#
			# list
			#
			for i, obj2 in enumerate(obj):
				if isinstance(obj2, list):
					#
					# nested list
					#
					for ii,jj in enumerate(obj2):
						jj = jj.strip()
					newobj[i] = jj
				elif isinstance(obj2, str):
					#
					# simple list
					#
					obj2 = obj2.replace(" ","")
					newobj[i] = obj2
				elif isinstance(obj2, int):
					#
					# integer
					#
					newobj[i] = obj2
				else:
					newobj[i] = obj2
		else: # error
			print "remove_space: input str or list"

		return newobj

def get_reac_and_prod(reactionfile):
	import numpy as np
	import os,sys
	#
	# form reactant and product information
	#
	(reac, rxn, prod) = read_reactionfile(reactionfile)

	rxn_num = len(rxn)

	r_ads  = [ [] for i in range(rxn_num) ]
	r_site = [ [] for i in range(rxn_num) ]
	r_coef = [ [] for i in range(rxn_num) ]

	p_ads  = [ [] for i in range(rxn_num) ]
	p_site = [ [] for i in range(rxn_num) ]
	p_coef = [ [] for i in range(rxn_num) ]


	for irxn, jrnx in enumerate(rxn):
		ireac = reac[irxn];     iprod = prod[irxn]
		ireac_num = len(ireac); iprod_num = len(iprod)
		#
		# reactant
		#
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
		#
		# product
		#
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

	return (r_ads, r_site, r_coef,  p_ads, p_site, p_coef)

def get_number_of_reaction(reactionfile):
	import numpy as np
	import os,sys
	#
	# return number of elementary reactions
	#
	(reac, rxn, prod) = read_reactionfile(reactionfile)
	rxn_num = len(rxn)
	return rxn_num


def get_preexponential(reactionfile):
	#
	# ! not needed for MATLAB use
	#
	from ase import Atoms, Atom
	from ase.collections import methane
	import numpy as np
	import os,sys
	#
	# calculate pre-exponential factor
	#
	(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

	rxn_num = get_number_of_reaction(reactionfile)

	Afor = np.array(rxn_num, dtype="f")
	Arev = np.array(rxn_num, dtype="f")

	mass_reac = np.array(rxn_num*[range(len(r_ads[0]))], dtype="f")
	mass_prod = np.array(rxn_num*[range(len(p_ads[0]))], dtype="f")

	for irxn in range(rxn_num):
		print "---------", irxn, "------------"
		#
		# reactants
		#
		for imol, mol in enumerate(r_ads[irxn]):
			tmp  = methane[mol]
			mass = sum(tmp.get_masses())
			mass_reac[irxn,imol] = mass
		#
		# products
		#
		for imol, mol in enumerate(p_ads[irxn]):
			tmp  = methane[mol]
			mass = sum(tmp.get_masses())
			mass_prod[irxn,imol] = mass

		Afor[irxn] = 1.0
		Arev[irxn] = 1.0

	return Afor,Arev

def get_rateconstant(reactionfile,Afor,Arev,Efor,Erev,T):
	#
	# ! not needed for MATLAB use
	#
	from ase import Atoms, Atom
	from ase.collections import methane
	import numpy as np
	import os,sys
	#
	# calculate rate constant
	#
	(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

	rxn_num = get_number_of_reaction(reactionfile)

	kfor = np.array(rxn_num, dtype="f")
	krev = np.array(rxn_num, dtype="f")

	RT = 8.314*T / 1000.0 # in kJ/mol

	for irxn in range(rxn_num):
		kfor[irxn] = Afor[irxn]*exp(-Efor[irxn] / RT)
		krev[irxn] = Arev[irxn]*exp(-Erev[irxn] / RT)

	return kfor,krev


def read_speciesfile(speciesfile):

	f = open(speciesfile)

	species = f.read()

	species = species.replace('[','')
	species = species.replace(']','')
	species = species.replace(' ','')
	species = species.replace('\'','')
	species = species.split(",")

	return species

def remove_parentheses(file):
	import os
	#
	# remove parentheses -- maybe for MATLAB use
	#
	tmpfile = "ttt.txt"
	os.system('cat %s | sed "s/\[//g" > %s' % (file, tmpfile) )
	os.system('cat %s | sed "s/\]//g" > %s' % (tmpfile, file) )
	os.system('rm %s' % tmpfile)


def get_species_num(species):
	from reaction_tools import read_speciesfile

	speciesfile = "species.txt"
	lst = read_speciesfile(speciesfile)

	return lst.index(species)


def get_adsorption_sites(infile):
	from reaction_tools import remove_space

	f = open(infile,"r")

	lines = f.readlines()

	mol  = range(len(lines))
	site = range(len(lines))

	for i,line in enumerate(lines):
		aaa,bbb = line.replace("\n","").split(":")
		mol[i]  = remove_space(aaa)
		bbb     = remove_space(bbb)
		site[i] = bbb.split(",")

	return mol,site

