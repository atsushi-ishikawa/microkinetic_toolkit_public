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

label_rxn = False

os.system('grep -v "^#" %s > reaction2.txt' % inp )
os.system('grep -v "^\s*$"   reaction2.txt > reaction3.txt')

numlines = sum(1 for line in open("reaction3.txt"))

f = open("reaction3.txt","r")
lines = f.readlines()
os.system('rm reaction2.txt reaction3.txt')

reac   = range(numlines)
rxn    = range(numlines)
prod   = range(numlines)
values = range(numlines)

edge_scale = 1.5
thre = 1.0

for i,line in enumerate(lines):
	if ':' in line:
		comp,value = line.split(':')
		value = value.replace('\n','')
		value = float(value)
		value = value if value > thre else 20.0
		value = log10(float(value))
	else:
	 	comp  = line
	 	value = 1.0

	comp = comp.replace('\n','').replace('>','').replace(' ','').split('--')
	reac_tmp = comp[0]
	rxn_tmp  = comp[1]
	prod_tmp = comp[2]

	reac[i]   = reac_tmp.split("+")
	#rxn[i]    = reac[i][0] + "_" + rxn_tmp
	rxn[i]    = 'rxn' + str(i+1)
	prod[i]   = prod_tmp.split("+")
	values[i] = value*edge_scale

c_siz = 200; c_col = "blue"
r_siz = 20;  r_col = "black"

if coverage:
	cov_dict = {}
	fcov = open(cov_file,"r")
	lines = fcov.readlines()
	for iline,line in enumerate(lines):
		cov = line.replace(' ','').replace('\n','')
		cov_dict[iline] = float(cov)

nodeA = 200.0
nodeB = 20.0

G = nx.Graph()
for i,j in enumerate(rxn):
	G.add_node(rxn[i], size=r_siz, color=r_col, typ='rxn')
 	for ireac,j1 in enumerate(reac[i]):
		if coverage:
			mol  = reac[i][ireac]
			spe  = get_species_num(mol)
			size = cov_dict[spe]
			size = nodeA*(nodeB + log10(size))
			size = int(size)
		else:
			size = c_siz

 		G.add_node(reac[i][ireac], size=size, color=c_col, typ='comp')
 		G.add_edge(reac[i][ireac], rxn[i], weight=values[i])
 
 	for iprod,j2 in enumerate(prod[i]):
		if coverage:
			mol  = reac[i][ireac]
			spe  = get_species_num(mol)
			size = cov_dict[spe]
			size = nodeA*(nodeB + log10(size))
			size = int(size)
		else:
			size = c_siz

 		G.add_node(prod[i][iprod], size=size, color=c_col, typ='comp')
 		G.add_edge(prod[i][iprod], rxn[i], weight=values[i])

#
# drawing 
#
siz = nx.get_node_attributes(G,'size')
col = nx.get_node_attributes(G,'color')

pos = nx.nx_pydot.graphviz_layout(G, prog='fdp') # prog='neato' is also a good choice
nx.draw_networkx_nodes(G, pos, nodelist=siz.keys(), node_size=siz.values(), node_color=col.values(), alpha=0.8)
edges = G.edges()
weights = [G[u][v]['weight'] for u,v in edges]

nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.6, width=weights)

# compound labels
Gcomp = nx.Graph()
for n,typ in G.nodes.data('typ'):
	if typ == 'comp':
		Gcomp.add_node(n)

nx.draw_networkx_labels(Gcomp, pos, font_size=14, font_family='Gill Sans MT')

# reaction labels
if label_rxn:
	Grxn = nx.Graph()
	for n,typ in G.nodes.data('typ'):
		if typ == 'rxn':
			Grxn.add_node(n)

	nx.draw_networkx_labels(Grxn, pos, font_size=12, font_family='Gill Sans MT')

plt.xticks([])
plt.yticks([])
plt.show()

# plt.figure(figsize=(16,10))
# plt.savefig("oxidative_coupling.eps", format="eps")

nx.write_gexf(G,'test.gexf')

