import os, sys
import matplotlib.pyplot as plt
import networkx as nx
from math import log10
from reaction_tools import *

argvs = sys.argv

if len(argvs)==3: # read coverage
	coverage = True
	cov_file = argvs[2]
else:
	coverage = False

inp = argvs[1]

os.system('grep -v "^#" %s > reaction2.txt' % inp )
os.system('grep -v "^\s*$"   reaction2.txt > reaction3.txt')

label_rxn = True

numlines = sum(1 for line in open("reaction3.txt"))

f = open("reaction3.txt","r")
lines = f.readlines()
os.system('rm reaction2.txt reaction3.txt')

reac   = range(numlines)
rxn    = range(numlines)
prod   = range(numlines)
values = range(numlines)

edge_scale = 1.0

for i,line in enumerate(lines):
	try:
		comp,value = line.split(':')
		value = value.replace('\n','')
		value = log10(float(value))
	except:
		comp  = line
		value = 1.0

	comp = comp.replace('\n','').replace('>','').replace(' ','').split('--')
	reac_tmp = comp[0]
	rxn_tmp  = comp[1]
	prod_tmp = comp[2]

	reac[i]   = reac_tmp.split("+")
	rxn[i]    = reac[i][0] + "_" + rxn_tmp
	prod[i]   = prod_tmp.split("+")
	values[i] = value*edge_scale

# G = nx.DiGraph()
G = nx.Graph()

c_siz = 200; c_col = "black"
r_siz = 20;  r_col = "black"

cov_scale = 400
if coverage:
	cov_dict = {}
	fcov = open(cov_file,"r")
	lines = fcov.readlines()
	for iline,line in enumerate(lines):
		cov = line.replace(' ','').replace('\n','')
		cov_dict[iline] = float(cov)*cov_scale

for i,j in enumerate(rxn):
	G.add_node(rxn[i], size=r_siz, color=r_col, typ="rxn")
 	for ireac,j1 in enumerate(reac[i]):
		if coverage:
			mol = reac[i][ireac]
			spe = get_species_num(mol)
			size = cov_dict[spe]
		else:
			size = c_siz
 		G.add_node(reac[i][ireac], size=size, color=c_col, typ="comp")
 		G.add_edge(reac[i][ireac], rxn[i], weight=values[i])
 
 	for iprod,j2 in enumerate(prod[i]):
		if coverage:
			mol = reac[i][ireac]
			spe = get_species_num(mol)
			size = cov_dict[spe]
		else:
			size = c_siz
 		G.add_node(prod[i][iprod], size=size, color=c_col, typ="comp")
 		G.add_edge(prod[i][iprod], rxn[i], weight=values[i])

#Greac = (n for n, d in G if G.nodes_iter(data=True)  if d["typ"]=="rxn")
#Gcomp = (n for n, d in G if G.nodes_iter(data=True)  if d["typ"]=="comp")

Greac = nx.Graph()
Gcomp = nx.Graph()

Greac = (n for n, d in G.nodes if G.nodes[d]['typ']=='rxn')
Gcomp = (n for n, d in G.nodes if G.nodes[d]['typ']=='comp')

#Greac.add_nodes_from([n for n, d in G.nodes_iter(data=True) if d["typ"]=="rxn"])
#Gcomp.add_nodes_from([n for n, d in G.nodes_iter(data=True) if d["typ"]=="comp"])

edges = G.edges()

#plt.figure(num=None, figsize=(20,20), dpi=350)
plt.figure(figsize=(16,10))

#pos = nx.spring_layout(G, k=0.4)
pos = nx.nx_pydot.graphviz_layout(G, prog="fdp") # prog="neato" is also good

siz = nx.get_node_attributes(G,"size")
col = nx.get_node_attributes(G,"color")

nx.draw_networkx_nodes(G, pos, nodelist=siz.keys(), node_size=siz.values(), node_color=col.values(),alpha=0.6)

#nx.draw_networkx_edges(G, pos, edge_color="gray",  alpha=1.0, width=0.6)
weights = [G[u][v]['weight'] for u,v in edges]
nx.draw_networkx_edges(G, pos, edge_color="gray",  alpha=1.0, width=weights)

if label_rxn:
	nx.draw_networkx_labels(G, pos, font_size=14, font_family="Gill Sans MT")
else:
 	nx.draw_networkx_labels(Gcomp, pos, font_size=14, font_family="Gill Sans MT")

plt.xticks([])
plt.yticks([])
plt.savefig("oxidative_coupling.eps", format="eps")
plt.show()

