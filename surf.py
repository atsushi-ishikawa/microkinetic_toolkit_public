from ase.build import fcc111, surface, sort, niggli_reduce, add_adsorbate
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read,write
from ase.visualize import view
from ase import Atoms, Atom
import os
import numpy as np
import collections
from reaction_tools import *

vacuum = 10.0
doping = False

if os.path.exists('surf.db'):
	os.system('rm surf.db')

nlayer = 3
nrelax = [1,2] # for each element

#cif_file = "ga2o3.cif"
#cif_file = "FePO4_P31_2_1.cif"
#cif_file = "LaCuO3.cif"
#cif_file = "SrCoO3.cif"
#cif_file = "LaMnO3.cif"
#cif_file = "wo3.cif"
#cif_file = "pd.cif"
#cif_file = "cao.cif"
cif_file = "La2O3.cif"
#cif_file = "Ce2W3O12.cif"

#lattice = "fcc"
facet   = "100"
lattice = "hcp"
#facet   = "001"
#lattice = "sp15"
#facet = "010"

indices = []
for c in facet:
	indices.append(int(c))
bulk = read(cif_file)
surf = surface(lattice=bulk, indices=indices, layers=nlayer, vacuum=vacuum) # step: (310) is good. nlayer=7, [1,2,1] might be good.

if cif_file == "La2O3.cif":
	surf.rotate(180,'y', rotate_cell=False) # La2O3
	surf.wrap()
	surf.center(axis=2) # La2O3, only z-axis

surf = surf*[3,2,1]
surf = sort(surf)
surf, zcount = sort_atoms_by_z(surf)

formula = surf.get_chemical_formula()

#offset_fac = 3.0/4.0
#offset_fac = 1.5 # for Pd110*[2,2,1]
offset_fac = (2.1,1.5)
offset_fac = np.array(offset_fac)
#offset_fac = 0.67
#
# doping e.g.) Mg by Li
#
if doping:
	symbols  =  np.array(surf.get_chemical_symbols())
	rep_atom = 48
	symbols[rep_atom] = 'Li'
	surf.set_chemical_symbols(symbols)

surf.translate([0,0,-vacuum+1])
#
# information for JSON file
#
pos = {	'lattice'   : lattice, 
		'facet'     : facet  ,
		'formula'   : formula,
		'offset_fac'  : offset_fac
      }
#
# relaxation issues
#
#natoms = len(surf.get_atomic_numbers())
symbols  = surf.get_chemical_symbols()
symbols  = sorted(symbols, key=symbols.index)
natoms   = collections.Counter(symbols)
elements = sorted(set(symbols), key=symbols.index)

last_num = []
tmp = 0
for i in elements:
	tmp += natoms[i]
	last_num.append(tmp)

#per_layer = natoms // nlayer
#print(per_layer)
#
# set tags: 2 = fixed layers, 1 = relaxing layers, 0 = adsorbates
#
# constraint will be set in reaction_energy.py
#
tag = np.full(sum(natoms.values()), 2, dtype=int)

if all([x==0 for x in nrelax]):
	print("hello")
#	for i in range(natoms):
#		tag[i] = 1
else:
	#
	# loop over each element group
	#
	for ielem, elem in enumerate(elements):
		now = last_num[ielem]
		nlayer_of_elem = len(zcount[ielem])
		for i in range(nlayer_of_elem, nlayer_of_elem-nrelax[ielem], -1):
			for j in range(now,now-zcount[ielem][i-1],-1):
				tag[j-1] = 1
				now -= 1

surf.set_tags(tag)

db = connect("surf.db")
db.write(surf, data=pos)
#
# adsorbate check
#
check_adsorbate = False
if check_adsorbate:
	mol = Atoms("H",[(0,0,0)])
	offset = (0.16, 0.33) # x1y1
	offset = np.array(offset)
	add_adsorbate(surf, mol, 0.6, position=(0,0), offset=offset*offset_fac)
	print(offset*offset_fac)

write("POSCAR",surf)

view(surf)

