import numpy as np
import os, sys, re, json, math
from reaction_tools import *
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.collections import methane
from ase.optimize import BFGS, MDMin, FIRE
from ase.vibrations import Vibrations, Infrared
from ase.db import connect
from ase.io import read, write
from ase.build import add_adsorbate, sort
from ase.visualize import view
from ase.neb import NEB, NEBTools
# -------------------------------------------------
# calculate reaction energy.
# molecule's data should be stored in "methane.json"
# -------------------------------------------------
# settings
#
argvs = sys.argv
reactionfile = argvs[1]

calculator = "vasp"
calculator = calculator.lower()
#
# temprary database to avoid overlapping calculations
#
dbfile = 'tmp.db'
dbfile = os.path.join(os.getcwd(), dbfile)
tmpdb  = connect(dbfile)
#
# if surface present, provide surface file
# in ase.db form
#
surface = True

if surface:
	dbfile = 'surf.db'
	dbfile = os.path.join(os.getcwd(), dbfile)
	db     = connect(dbfile)
	surf       = db.get_atoms(id=1)
	lattice    = db.get(id=1).data.lattice
	facet      = db.get(id=1).data.facet
	surf_name  = db.get(id=1).data.formula
	offset_fac = db.get(id=1).data.offset_fac

	# load site information
	f = open('site_info.json','r')
	site_info = json.load(f)
	f.close()

# fix atoms
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 2])
surf.set_constraint(c)

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

rxn_num = get_number_of_reaction(reactionfile)

maxoptsteps = 200
ads_height0 = 1.6
ads_pos0 = (0.0, 0.0)

ZPE = [False, False]
IR  = [False, False] # whether to do IR...[Reac, Prod]

# transition state
TS = True

# single point
SP = False

if TS:
	CI = False # whether to do CI-NEB

nimages = 8

if TS:
	vtst = "/home/a_ishi/vasp/vtstscripts/vtstscripts-935/" # whisky
	if not os.path.exists(vtst):
		vtst = "/home/usr6/m70286a/vasp/vtstscripts/vtstscripts-935/" # kyushu

# whether to do single point after optimization
# at different computational level

## --- Gaussian ---
if "gau" in calculator:
	method = "b3lyp"
	basis  = "6-31G(d)" # do not use aesterisk for polarization func
	nprocs = 12
	mem    = "8GB"
	if SP:
		method_sp = "ccsd(t)"
	basis_name = re.sub("\(", "", basis)
	basis_name = re.sub("\)", "", basis_name)
	basis_name = re.sub(",",  "", basis_name)
	label = method + "_" + basis_name

## --- VASP ---
elif "vasp" in calculator:
	xc          = "beef-vdw"
	ivdw        = 0
	# GGA list
	#  GGAs: pw91, pbe, pbesol, revpbe, rpbe, am05
	#  meta-GGAs: tpss, revtpss, m06l, ms0, ms1, scan, scan-rvv10
	#    --> gga and pp (to be override) are set automatically
	#  vdw-DFs: vdw-df, optpbe-vdw, optb88-vdw, optb86b-vdw, vdw-df2, beef-vdw
	#    --> luse_vdw and others are set automatically
	prec        = "normal"
	encut       = 400.0 # 213.0 or 400.0 or 500.0
	potim       = 0.10
	ibrion      = 1
	nfree       = 10
	nsw         = 200
	nsw_neb     = 1
	nsw_dimer   = 1
	nelmin      = 5
	nelm        = 100 # default:40
	ediff       = 1.0e-6
	ediffg      = -0.05
	kpts_surf   = [5, 5, 1]
	ismear_surf = 1
	sigma_surf  = 0.10
	vacuum      = 10.0 # for gas-phase molecules. surface vacuum is set by surf.py
	setups      = None
	ialgo       = 48 # normal=38, veryfast=48
	npar        = 6
	nsim        = npar
	lwave       = False
	lcharg      = False
	ispin       = 2
	#setups = {"O" : "_h"}

	# set lmaxmix
	lmaxmix = 2
	first_TM = ["Fe","Ni"]
	for elem in first_TM:
		if elem in surf.get_chemical_symbols():
			lmaxmix = 4

	if xc=='pw91':
		pp = 'potpaw_GGA'
	else: # mostly PBE is OK
		pp = 'potpaw_PBE.54'

	# dipoole correction
	dipole = True
	if dipole:
		ldipol = True
		idipol = 3 # xyz-direction

	method = xc
	basis  = ""
	label  = method

	# DFT+U
	DFTU = False
	if DFTU:
		ldau = "true"
		ldautype = 2
		ldau_luj = { 'La':{'L':3, 'U':3.5, 'J':0.0}, 'O':{'L':-1, 'U':0.0, 'J':0.0} }
		ialgo = 38

	# charge
	neutral = True
	if neutral:
		charge = 0
	else:
		charge = 1

	if SP:
		xc_sp = "hse06"
		if ("hse" or "b3lyp") in xc_sp:
			ialgo_sp = 58
		ivdw_sp = 0

## --- EMT --- -> nothing to set
elif "emt" in calculator:
	method = ""
	basis  = ""
	label  = ""

if True in ZPE:
	label = label + "_ZPE"
if SP:
	label = label + "_SP"

barrierfile  = reactionfile.split(".")[0] + "_Ea_" + label + ".txt"
deltaEfile = "deltaE.txt"
fbarrier = open(barrierfile, 'w')
fbarrier.close()
fdeltaE = open(deltaEfile, 'w')
fdeltaE.close()

print "calculator:" + calculator + " method: " + method + " basis: " + basis

for irxn in range(rxn_num):
	fbarrier = open(barrierfile, 'a')
	fdeltaE  = open(deltaEfile,  'a')
	print "--- calculating elementary reaction No. ", irxn, "---"

	reac_en = np.array(range(len(r_ads[irxn])),dtype="f")
	prod_en = np.array(range(len(p_ads[irxn])),dtype="f")
	#
	# reactants
	#
	for imols, mols in enumerate(r_ads[irxn]):
		surf_tmp = surf.copy()

		for imol, mol in enumerate(mols):
			ads_height = ads_height0
			print "----- reactant: molecule No.", imol, " is ", mol, "-----"
			config = "normal"

			if 'surf' in mol:
				# surface itself
				mol,neutral,charge = read_charge(mol)
				tmp = surf
			elif 'def' in mol:
				# defected surface
				mol,neutral,charge = read_charge(mol)
				tmp = 'def'
			else:
				# some adsorbate is present
				mol,neutral,charge = read_charge(mol)

				# flip or rotate
				if "-SIDEy" in mol:
					# if "-SIDE" only, maybe this direction
					mol = mol.replace("-SIDEy","")
					tmp = methane[mol]
					tmp.rotate(90,'y')
					config = "sidey"
				elif "-SIDEx" in mol:
					mol = mol.replace("-SIDEx","")
					tmp = methane[mol]
					tmp.rotate(90,'x')
					config = "sidex"
				elif "-FLIP" in mol:
					mol = mol.replace("-FLIP","")
					tmp = methane[mol]
					tmp.rotate(180,'y')
					config = "flip"
				elif "-TILT" in mol:
					mol = mol.replace("-TILT","")
					tmp = methane[mol]
					tmp.rotate(45,'y')
					config = "tilt"
				elif "-HIGH" in mol:
					mol = mol.replace("-HIGH","")
					tmp = methane[mol]
					ads_height += 1.5
					config = "high"
				elif "-AIR" in mol:
					mol = mol.replace("-AIR","")
					tmp = methane[mol]
					ads_height += 4.0
					config = "air"
				else:
					tmp = methane[mol]

			site = r_site[irxn][imols][imol]

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			if site == 'gas':
				if 'surf' in mols:
					mol_type = 'surf'
				else:
					mol_type = 'gaseous'
			else:
				mol_type = 'adsorbed'
			
			if mol_type=='adsorbed':
				offset = site_info[lattice][facet][site][site_pos]
				if len(offset)==3:
					shift = offset[2] * surf.get_cell()[2][2]
					ads_height += shift
					offset = offset[0:2]

				offset = np.array(offset)
				#
				# wrap atoms to prevent adsorbate being on different cell
				#
				surf_tmp.translate([0,0,2])
				surf_tmp.wrap(pbc=[0,0,1])
				surf_tmp.translate([0,0,-2])
				print("lattice:{0}, facet:{1}, site:{2}, site_pos:{3}, config:{4}".format(lattice, facet, site, site_pos, config))
				#
				if tmp == 'def':
					defect = find_closest_atom(surf_tmp, offset=offset*offset_fac)
					del surf_tmp[len(surf_tmp.get_atomic_numbers())-1]
					del surf_tmp[defect] # defect
					tmp = surf_tmp
				else:
					#
					# shift adsorbate molecule
					#
					xcom  = tmp.get_center_of_mass()[0]
					ycom  = tmp.get_center_of_mass()[1]
					shift = (xcom - tmp.positions[0,0], ycom - tmp.positions[0,1])
					z_shift = tmp.positions[:,2].min()

					if tmp.get_chemical_formula()  == 'H': # special attention to H
						ads_height = 1.2

					ads_height -= z_shift
					ads_pos = (ads_pos0[0]-shift[0], ads_pos0[1]-shift[1])
					add_adsorbate(surf_tmp, tmp, ads_height, position=ads_pos, offset=offset*offset_fac)
					tmp = surf_tmp
		del surf_tmp
		#view(tmp); quit()
		#
		# end adsorbing molecule
		#

		#
		# Identification done. Look for temporary database 
		# for identical system.
		#
		formula = tmp.get_chemical_formula()
		try:
			past = tmpdb.get(formula=formula)
		except:
			print "first time"
			first_time = True
		else:
			if site == past.data.site:
 				if site_pos == past.data.site_pos:
 					if config == past.data.config:
						if len(mols) == 1:
 							print "already calculated"
 							tmp = tmpdb.get_atoms(id=past.id)
 							first_time = False

		magmom = tmp.get_initial_magnetic_moments()
		natom  = len(tmp.get_atomic_numbers())
		coef   = r_coef[irxn][imols]

		# spin switch
		if int( math.ceil(sum(magmom)) )!=0:
			print "has unpaired electron"
			ispin = 2
		else:
			ispin = 1

		if natom == 1:
			ZPE[0] = False; IR[0] = False
		#
		# set label
		#
		r_label = label + "_rxn" + str(irxn) + "_" + "-".join(mols) + "_" + site
		if mol_type=='adsorbed':
			r_label = r_label + "_" + surf_name
		r_traj  = r_label + "reac.traj"
		#
		# branch compurational setting by gas or not
		# 
		if site == 'gas' and not 'surf' in mols:
			#
 			# gas-phase molecule
			#
			if 'vasp' in calculator:
				cell = np.array([1, 1, 1])
				cell = vacuum*cell
				tmp.set_cell(cell)
				tmp.center()
				ismear = 0 # gaussian smearing
				sigma  = 0.05
				kpts = [1,1,1]
				gas_mol = True
		else: # surface
			ismear  = ismear_surf
			sigma   = sigma_surf
			kpts    = kpts_surf
			gas_mol = False
		#
		# set calculator
		#
		if "gau" in calculator:
			tmp.calc = Gaussian(label=r_label, method=method, basis=basis, scf="xqc", nprocs=nprocs, mem=mem)
			opt = BFGS(tmp, trajectory=r_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)
			if SP:
				r_label = r_label + "_sp"
				tmp.calc = Gaussian(label=r_label, method=method_sp, basis=basis, force=None, nprocs=nprocs, mem=mem)
		elif "vasp" in calculator:
			if gas_mol:
				#
				# gas-phase molecules
				#
				if neutral:
		 			tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol, idipol=idipol )
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
		 			tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
									nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol )
			else:
				#
				# surface
				#
				if neutral:
					if DFTU:
		 				tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj, pp=pp, ldipol=ldipol, idipol=idipol ) # DFT+U
					else:
			 			tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol, 
										idipol=idipol, nfree=nfree ) # normal
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
					if DFTU:
		 				tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, nelect=nelect, lmono="true", 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj, pp=pp, ldipol=ldipol, idipol=idipol ) # charge + DFT+U
					else:
		 				tmp.calc = Vasp(label=r_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol ) # charge

		elif "emt" in calculator:
			tmp.calc = EMT()
			opt = BFGS(tmp, trajectory=r_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)

		en = tmp.get_potential_energy()

		# single point
		if SP:
			if mol_type=='gaseous':
				ismear_sp = 0
				kpts_sp   = [1, 1, 1]
			else:
				#ismear_sp = -5
				ismear_sp = ismear
				#kpts_sp   = [5, 5, 1]
				kpts_sp = kpts

			tmp.calc = Vasp(label=r_label, prec=prec, xc=xc_sp, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw_sp, npar=npar, nsim=nsim,
							encut=encut, ismear=ismear_sp, istart=0, setups=setups, sigma=sigma, ialgo=ialgo_sp, lwave=lwave, lcharg=lcharg,
							ibrion=-1, potim=potim, nsw=0, ediff=ediff, ediffg=ediffg, kpts=kpts_sp, pp=pp, ldipol=ldipol, idipol=idipol) # normal
			en = tmp.get_potential_energy()

		if "vasp" in calculator:
			xmlfile = "vasprun_" + r_label + ".xml"
			contcar = "CONTCAR_" + r_label
			reac_contcar = contcar
			os.system('cp vasprun.xml %s >& /dev/null' % xmlfile)
			os.system('cp CONTCAR %s >& /dev/null' % contcar)
			os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHGCAR CHG WAVECAR REPORT >& /dev/null')

		if mol_type!='surf' and (ZPE[0] or IR[0]): # surf--nothing to do with vibration
			# fix atoms for vibrations
			if mol_type=='adsorbed':
				c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1 or atom.tag == 2])
				tmp.set_constraint(c)
			if ZPE[0]:
				# do ZPE calculation
				zpedir = r_label + "_ZPE"
				if not os.path.isdir(zpedir):
					os.makedirs(zpedir)
				vib = Vibrations(tmp, delta=0.01, name=zpedir+"/"+"ZPE") # delta = 0.01 is default
				vib.run()
				hnu = vib.get_energies()
				zpe = vib.get_zero_point_energy()
				en = en + zpe
				os.system("rm -r zpedir")
			if IR[0]:
				# setting for IR calculation
				irdir = r_label + "_IR"
				if not os.path.isdir(irdir):
					os.makedirs(irdir)

				tmp.calc = Vasp(label=irdir, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								ediff=ediff, kpts=kpts, isym=0, lmono=True, ldipol=ldipol, idipol=idipol, dipol=tmp.get_center_of_mass(scaled=True), pp=pp )
				vib = Infrared(tmp, delta=0.01, name=irdir+"/"+"IR")
				vib.run()
				# vib.write_spectra(out=r_label+"_IR.dat",start=1000,end=4000, width=10, normalize=True)
				vib.write_spectra(out=irdir+"_spectra.dat", start=1000, end=4000)
				vib.write_mode()
				vib.write_jmol()
				# include ZPE
				zpe = vib.get_zero_point_energy()
				en = en + zpe
				os.system("rm *.pckl")

			reac_en[imols] = en
		else:
			reac_en[imols] = en

		reac_en[imols] = coef * reac_en[imols]

		# recording to database
		if(first_time):
			tmpdb.write(tmp, data={'site':site, 'site_pos':site_pos, 'config':config})

	#
	# products
	#
	for imols, mols in enumerate(p_ads[irxn]):
		surf_tmp = surf.copy()

		for imol, mol in enumerate(mols):
			ads_height = ads_height0
			print "----- product: molecule No.", imol, " is ", mol, "-----"
			config = "normal"

			if 'surf' in mol:
				# surface itself
				mol,neutral,charge = read_charge(mol)
				tmp = surf
			elif 'def' in mol:
				# defected surface
				mol,neutral,charge = read_charge(mol)
				tmp = 'def'
			else:
				# some adsorbate is present
				mol,neutral,charge = read_charge(mol)

				# flip or rotate
				if "-SIDEy" in mol:
					mol = mol.replace("-SIDEy","")
					tmp = methane[mol]
					tmp.rotate(90,'y')
					config = "sidey"
				elif "-SIDEx" in mol:
					mol = mol.replace("-SIDEx","")
					tmp = methane[mol]
					tmp.rotate(90,'x')
					config = "sidex"
				elif "-FLIP" in mol:
					mol = mol.replace("-FLIP","")
					tmp = methane[mol]
					tmp.rotate(180,'y')
					config = "flip"
				elif "-TILT" in mol:
					mol = mol.replace("-TILT","")
					tmp = methane[mol]
					tmp.rotate(45,'y')
					config = "tilt"
				elif "-HIGH" in mol:
					mol = mol.replace("-HIGH","")
					tmp = methane[mol]
					ads_height += 1.5
					config = "high"
				elif "-AIR" in mol:
					mol = mol.replace("-AIR","")
					tmp = methane[mol]
					ads_height += 4.0
					config = "air"
				else:
					tmp = methane[mol]

			site = p_site[irxn][imols][imol]

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			if site == 'gas':
				if 'surf' in mols:
					mol_type = 'surf'
				else:
					mol_type = 'gaseous'
			else:
				mol_type = 'adsorbed'

			if mol_type=='adsorbed':
				offset = site_info[lattice][facet][site][site_pos]
				if len(offset)==3:
					shift = offset[2] * surf.get_cell()[2][2]
					ads_height += shift
					offset = offset[0:2]

				offset = np.array(offset)
				#
				# wrap atoms to prevent adsorbate being on different cell
				#
				surf_tmp.translate([0,0,2])
				surf_tmp.wrap(pbc=[0,0,1])
				surf_tmp.translate([0,0,-2])
				print("lattice:{0}, facet:{1}, site:{2}, site_pos:{3}, config:{4}".format(lattice, facet, site, site_pos, config))
				#
				if tmp == 'def':
					defect = find_closest_atom(surf_tmp, offset=offset*offset_fac)
					del surf_tmp[len(surf_tmp.get_atomic_numbers())-1]
					del surf_tmp[defect] # defect
					tmp = surf_tmp
				else:
					#
					# shift adsorbate molecule
					#
					xcom  = tmp.get_center_of_mass()[0]
					ycom  = tmp.get_center_of_mass()[1]
					shift = (xcom - tmp.positions[0,0], ycom - tmp.positions[0,1])
					z_shift = tmp.positions[:,2].min()

					if tmp.get_chemical_formula()  == 'H': # special attention to H
						ads_height = 1.2

					ads_height -= z_shift
					ads_pos = (ads_pos0[0]-shift[0], ads_pos0[1]-shift[1])
					add_adsorbate(surf_tmp, tmp, ads_height, position=ads_pos, offset=offset*offset_fac)
					tmp = surf_tmp
		del surf_tmp
		#
		# end adsorbing molecule
		#

		#
		# Identification done. Look for temporary database 
		# for identical system.
		#
		formula = tmp.get_chemical_formula()
		try:
			past = tmpdb.get(formula=formula)
		except:
			print "first time"
			first_time = True
		else:
			if site == past.data.site:
				if site_pos == past.data.site_pos:
					if config == past.data.config:
						if len(mols) == 1:
							print "already calculated"
							tmp = tmpdb.get_atoms(id=past.id)
							first_time = False

		magmom = tmp.get_initial_magnetic_moments()
		natom  = len(tmp.get_atomic_numbers())
		coef   = p_coef[irxn][imols]

		# spin switch
		if int( math.ceil(sum(magmom)) )!=0:
			print "has unpaired electron"
			ispin = 2
		else:
			ispin = 1

		if natom == 1:
			ZPE[1] = False; IR[1] = False
		#
		# set label
		#
		p_label = label + "_rxn" + str(irxn) + "_" + "-".join(mols) + "_" + site
		if mol_type=='adsorbed':
			p_label = p_label + "_" + surf_name
		p_traj  = p_label + "prod.traj"
		#
		# branch compurational setting by gas or not
		# 
		if site == 'gas' and not 'surf' in mols:
			#
 			# gas-phase molecule
			#
			if 'vasp' in calculator:
				cell = np.array([1, 1, 1])
				cell = vacuum*cell
				tmp.set_cell(cell)
				tmp.center()
				ismear = 0 # gaussian smearing
				sigma  = 0.05
				kpts = [1,1,1]
				gas_mol = True
		else: # surface
			ismear  = ismear_surf # Methfessel-Paxton
			sigma   = sigma_surf
			kpts    = kpts_surf
			gas_mol = False
		#
		# set calculator
		#
		if "gau" in calculator:
			tmp.calc = Gaussian(label=p_label, method=method, basis=basis, scf="xqc", nprocs=nprocs, mem=mem)
			opt = BFGS(tmp, trajectory=p_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)
			if SP:
				p_label = p_label + "_sp"
				tmp.calc = Gaussian(label=p_label, method=method_sp, basis=basis, force=None, nprocs=nprocs, mem=mem)
		elif "vasp" in calculator:
			if gas_mol:
				#
				# gas-phase molecules
				#
				if neutral:
		 			tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol, idipol=idipol )
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
		 			tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
									nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol )
			else:
				#
				# surface
				#
				if neutral:
					if DFTU:
		 				tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj, pp=pp, ldipol=ldipol, idipol=idipol ) # DFT+U
					else:
			 			tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol,
										idipol=idipol, nfree=nfree ) # normal
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
					if DFTU:
		 				tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, nelect=nelect, lmono="true", 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj, pp=pp, ldipol=ldipol, idipol=idipol ) # charge + DFT+U
					else:
		 				tmp.calc = Vasp(label=p_label, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
										ibrion=ibrion, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol ) # charge

		elif "emt" in calculator:
			tmp.calc = EMT()
			opt = BFGS(tmp, trajectory=p_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)

		en = tmp.get_potential_energy()

		# single point
		if SP:
			if mol_type=='gaseous':
				ismear_sp = 0
				kpts_sp   = [1, 1, 1]
			else:
				#ismear_sp = -5
				ismear_sp = ismear
				#kpts_sp   = [5, 5, 1]
				kpts_sp = kpts

			tmp.calc = Vasp(label=p_label, prec=prec, xc=xc_sp, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw_sp, npar=npar, nsim=nsim,
							encut=encut, ismear=ismear_sp, istart=0, setups=setups, sigma=sigma, ialgo=ialgo_sp, lwave=lwave, lcharg=lcharg,
							ibrion=-1, potim=potim, nsw=0, ediff=ediff, ediffg=ediffg, kpts=kpts_sp, pp=pp, ldipol=ldipol, idipol=idipol ) # normal
			en = tmp.get_potential_energy()

		if "vasp" in calculator:
			xmlfile = "vasprun_" + p_label + ".xml"
			contcar = "CONTCAR_" + p_label
			prod_contcar = contcar
			os.system('cp vasprun.xml %s >& /dev/null' % xmlfile)
			os.system('cp CONTCAR %s >& /dev/null' % contcar)
			os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHGCAR CHG WAVECAR REPORT >& /dev/null')

		if mol_type!='surf' and (ZPE[1] or IR[1]): # surf--nothing to do with vibration
			# fix atoms for vibrations
			if mol_type=='adsorbed':
				c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1 or atom.tag == 2])
				tmp.set_constraint(c)
			if ZPE[1]:
				# do ZPE calculation
				zpedir = p_label + "_ZPE"
				if not os.path.isdir(zpedir):
					os.makedirs(zpedir)
				vib = Vibrations(tmp, delta=0.01, name=zpedir+"/"+"ZPE") # delta = 0.01 is default
				vib.run()
				hnu = vib.get_energies()
				zpe = vib.get_zero_point_energy()
				en = en + zpe
				os.system("rm -r zpedir")
			if IR[1]:
				# setting for IR calculation
				irdir = p_label + "_IR"
				if not os.path.isdir(irdir):
					os.makedirs(irdir)

				tmp.calc = Vasp(label=irdir, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								ediff=ediff, kpts=kpts, isym=0, lmono=True, ldipol=ldipol, idipol=idipol, dipol=tmp.get_center_of_mass(scaled=True), pp=pp )
				vib = Infrared(tmp, delta=0.01, name=irdir+"/"+"IR")
				vib.run()
				# vib.write_spectra(out=r_label+"_IR.dat",start=1000,end=4000, width=10, normalize=True)
				vib.write_spectra(out=irdir+"_spectra.dat", start=1000, end=4000)
				vib.write_mode()
				vib.write_jmol()
				# include ZPE
				zpe = vib.get_zero_point_energy()
				en = en + zpe
				os.system("rm *.pckl")

			prod_en[imols] = en
		else:
			prod_en[imols] = en

		prod_en[imols] = coef * prod_en[imols]

		# recording to database
		if(first_time):
			tmpdb.write(tmp, data={'site':site, 'site_pos':site_pos, 'config':config})

		#
		# TS calc
		#
		if(TS):
			#
			# make different directory and go there
			#
			os.system('rm -rf tsdir') # remove old one

			if not os.path.exists("tsdir"):
				os.makedirs("tsdir")
			os.chdir("tsdir")

			# copy reactant and product POSCAR as POSCAR_reac and POSCAR_prod
			contcar1 = "../" + r_label + "/CONTCAR"
			contcar2 = "../" + p_label + "/CONTCAR"
			os.system('cp %s ./POSCAR_reac' % contcar1)
			os.system('cp %s ./POSCAR_prod' % contcar2)

			# Swith atomic index in order to make POSCAR_reac and POSCAR_prod close.
			# Different ordering cause bad NEB images.
			atom1 = read('POSCAR_reac')
			atom2 = read('POSCAR_prod')
			newatom1 = make_it_closer_by_exchange(atom1, atom2, thre=100.0) # atom1 is exchanged
			write('POSCAR_reac',newatom1) # do not sort because POTCAR does not follow
			write('POSCAR_prod',atom2)
			#
			# do "nebmake.pl"
			#
			nebmake = vtst + "nebmake.pl"

			os.system('%s POSCAR_reac POSCAR_prod %d >& /dev/null' % (nebmake,nimages))

			outcar1  = "../" + r_label + "/OUTCAR"
			outcar2  = "../" + p_label + "/OUTCAR"
			# copy reactant and product CONTCAR and OUTCAR to 00 and 0#IMAGE's POSCAR and OUTCAR
			os.system('cp %s 00'          % outcar1)
			os.system('cp %s 00/POSCAR'   % contcar1)
			os.system('cp %s %02d'        % (outcar2,  nimages+1))
			os.system('cp %s %02d/POSCAR' % (contcar2, nimages+1))

			# normal NEB
			if neutral:
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
			 					ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts, 
			 					images=nimages, spring=-5.0, lclimb=False, iopt=7, maxmove=0.10, pp=pp, ldipol=ldipol, idipol=idipol )
			else:
				nelect = get_number_of_valence_electrons(tmp)
				nelect = nelect - charge
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
			 					ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts, 
			 					images=nimages, spring=-5.0, lclimb=False, iopt=7, maxmove=0.10, nelect=nelect, lmono="true", 
								pp=pp, ldipol=ldipol, idipol=idipol )

			print "----------- doing NEB calculation with images=",nimages,"-----------"
			tmp.get_potential_energy()
			print "----------- normal NEB done -----------"
			neb_copy_contcar_to_poscar(nimages)

			# CI-NEB
			if CI:
				if neutral:
					tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts, 
									images=nimages, spring=-5.0, lclimb=True, iopt=7, maxmove=0.10, pp=pp, ldipol=ldipol, idipol=idipol )
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
					tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
									ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts, images=nimages,
									spring=-5.0, lclimb=True, iopt=7, maxmove=0.10, nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol )
				print "---------- doing CI-NEB calculation with images=",nimages,"-----------"
				tmp.get_potential_energy()
				print "----------- CI NEB done -----------"
				neb_copy_contcar_to_poscar(nimages)
 
			nebresults = vtst + "nebresults.pl"
			neb2dim    = vtst + "neb2dim.pl"
			os.system('%s >& /dev/null' % nebresults)
			os.system('%s >& /dev/null' % neb2dim)
			os.chdir("dim")

			# print "dimer made" ; quit()

			# dimer method
			if neutral:
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								nsw=nsw_dimer, ediff=ediff*0.1, ediffg=ediffg*0.5, kpts=kpts,
								iopt=2, maxmove=0.20, dfnmax=1.0, ichain=2, pp=pp, ldipol=ldipol, idipol=idipol )
			else:
				nelect = get_number_of_valence_electrons(tmp)
				nelect = nelect - charge
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
								encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								nsw=nsw_dimer, ediff=ediff*0.1, ediffg=ediffg*0.5, kpts=kpts,
								iopt=2, maxmove=0.20, dfnmax=1.0, ichain=2, nelect=nelect, lmono="true", pp=pp, ldipol=ldipol, idipol=idipol )

			print "----------- doing dimer method TS optimization -----------"
			TSene = tmp.get_potential_energy()
			print "----------- dimer done -----------"

			Ea = TSene -reac_en
			print "Ea = ",Ea

	deltaE = np.sum(prod_en) - np.sum(reac_en)
	print "deltaE=",deltaE
	#
	# writing reaction
	#
	string = ""
	for imol, mol in enumerate(r_ads[irxn]):
		mol_join = '-'.join(mol)
		string = string + "{0}_{1}".format(mol_join, '-'.join(r_site[irxn][imol]))
		if imol != len(r_ads[irxn])-1:
			string = string + " + "

	string = string + " --> "

	for imol, mol in enumerate(p_ads[irxn]):
		mol_join = '-'.join(mol)
		string = string + "{0}_{1}".format(mol_join, '-'.join(p_site[irxn][imol]))
		if imol != len(p_ads[irxn])-1:
			string = string + " + "

	fbarrier.write('{0:<70s}'.format(string))

	Eafor  =  deltaE
	Earev  = -deltaE
	Ea = [Eafor, Earev]

	fbarrier.write('{0:>14.8f} {1:>14.8f}\n'.format(Eafor, Earev))
	fbarrier.close()

	fdeltaE.write('{0:>14.8f} {1:>14.8f}\n'.format(Eafor, Earev))
	fdeltaE.close()
	#
	# loop over reaction
	#
remove_parentheses(barrierfile)
os.system("rm tmp.db >& /dev/null") # delte temporary database

