import numpy as np
import os, sys ,json
from reaction_tools import *
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.collections import methane
from ase.optimize import MDMin
from ase.vibrations import Vibrations
from ase.db import connect
from ase.io import read
from ase.build import add_adsorbate
from ase.visualize import view
# -------------------------------------------------
# calculate adsorption energy.
# molecule's data should be stored in "methane.json"
# -------------------------------------------------
argvs = sys.argv

ads_file   = argvs[1]
e_ads_file = ads_file.split(".")[0] + "_Ead" + ".txt"

calculator = "vasp" ; calculator = calculator.lower()

db_surf = connect('surf.db')
surf    = db_surf.get_atoms(id=1)
lattice = db_surf.get(id=1).data.lattice
facet   = db_surf.get(id=1).data.facet

# load site information
f = open('site_info.json','r')
site_info = json.load(f)
f.close()

# define json file
ads_json = 'e_ads.json'
db_ads = connect(ads_json)

fads = open(e_ads_file, "w")

mols,sites = get_adsorption_sites(ads_file)
mols_num = len(mols)

e_ads = np.array([0]*mols_num, dtype="f")

# parameters
ZPE = False
maxoptsteps = 50
ads_hight = 2.3

## --- VASP ---
if "vasp" in calculator:
	xc    = "pbe"
	prec  = "normal"
	encut = 400.0
	potim = 0.10
	nsw   = 100
	ediff = 1.0e-4
	ediffg = -0.03
	kpts = [3, 3, 1]
## --- EMT --- -> nothing to set

for imol, mol in enumerate(mols):
	print "----- reactant: molecule No.", imol, " is ", mol, "-----"

	tmp = methane[mol]

	for isite, site in enumerate(sites[imol]):
		surf_tmp = surf.copy()
		try:
			site_pos = site.split(".")[1]
		except:
			site_pos = "x1y1"

		print "lattice",lattice; print "facet", facet; print "site",site; print "site_pos",site_pos
		if site != 'gas':
			offset = site_info[lattice][facet][site][site_pos]
			add_adsorbate(surf_tmp, tmp, ads_hight, position=(0,0), offset=offset)
			tmp2 = surf_tmp
			del surf_tmp
		else:
			tmp2 = tmp
			tmp2.set_cell([10.0, 10.0, 10.0])

		magmom = tmp2.get_initial_magnetic_moments()
		natom  = len(tmp2.get_atomic_numbers())
		traj = mol + "_" + site + "_" + "ads.traj"

		if "vasp" in calculator:
			tmp2.calc = Vasp(prec=prec,xc=xc,ispin=2,encut=encut, ismear=0, istart=0,
						ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts )
		elif "emt" in calculator:
			tmp2.calc = EMT()
			opt = MDMin(tmp2, trajectory=traj)
			opt.run(fmax=0.05, steps=maxoptsteps)

		en  = tmp2.get_potential_energy()

		print "Eads = ", en

#		fads.write(mol + "_" + site + "   " + str(en) + "\n")
		fads.write('{0:<8}{1:<8}{2:<16.8f}\n'.format(mol,site,en))

		db_ads.write(tmp2, data={ mol + "-" + site: en})

		del tmp2

		if "vasp" in calculator:
			xmlfile = "vasprun_" + mol + "_" + site + ".xml"
			os.system("cp vasprun.xml %s",xmlfile)

		# loop over site
	# loop over molecule

fads.close()

