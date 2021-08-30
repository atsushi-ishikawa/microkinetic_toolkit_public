import numpy as np
import os, sys, re, json, math
import argparse
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
parser = argparse.ArgumentParser()
parser.add_argument("--reactionfile", required=True, help="file with elementary reactions")
parser.add_argument("--start", default=0, help="elementary reaction number to start")
parser.add_argument("--end",   default=0, help="elementary reaction number to end")
parser.add_argument("--calculator", default="emt", choices=["emt", "vasp"])
argvs = parser.parse_args()
reactionfile = argvs.reactionfile
rxnst = argvs.start
rxned = argvs.end
calculator = argvs.calculator.lower()

rxn_num = get_number_of_reaction(reactionfile)
if rxned == 0:  # not set...default
	rxned = rxn_num

clean = True  # whether to remove previous calculation
traj_dir = "traj_dir"
if clean:
	os.system("rm tmp.db")
	os.system("rm -rf %s" % traj_dir)
#
# temprary database to avoid overlapping calculations
#
tmpdbfile = "tmp.db"
tmpdbfile = os.path.join(os.getcwd(), tmpdbfile)
tmpdb = connect(tmpdbfile)
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
	f = open('site_info.json', 'r')
	site_info = json.load(f)
	f.close()

	# fix atoms
	c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 2])
	surf.set_constraint(c)

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

maxoptsteps = 200
ads_height0 = 1.8
ads_pos0 = (0.0, 0.0)

ZPE = {"reactant": False, "product": False}  # whether to include ZPE
IR  = {"reactant": False, "product": False}  # whether to do IR

# transition state
TS = False

# single point
SP = False

if TS:
	CI = False  # whether to do CI-NEB
	nimages = 8
	vtst = "/home/a_ishi/vasp/vtstscripts/vtstscripts-935/"  # whisky
	if not os.path.exists(vtst):
		vtst = "/home/usr6/m70286a/vasp/vtstscripts/vtstscripts-935/"  # kyushu

# whether to do single point after optimization
# at different computational level

## --- Gaussian ---
if calculator == "gaussian":
	method = "b3lyp"
	basis  = "6-31G(d)"  # do not use aesterisk for polarization func
	nprocs = 12
	mem    = "8GB"
	basis_name = re.sub("\(", "", basis)
	basis_name = re.sub("\)", "", basis_name)
	basis_name = re.sub(",",  "", basis_name)
	label = method + "_" + basis_name

	calc_gas = Gaussian(label=dir, method=method, basis=basis, scf="xqc", nprocs=nprocs, mem=mem)
	if SP:
		method_sp = "ccsd(t)"
		calc_gas_sp = Gaussian(label=dir, method=method_sp, basis=basis, scf="xqc", nprocs=nprocs, mem=mem)

## --- VASP ---
elif calculator == "vasp":
	# GGA list
	#  GGAs: pw91, pbe, pbesol, revpbe, rpbe, am05
	#  meta-GGAs: tpss, revtpss, m06l, ms0, ms1, scan, scan-rvv10
	#    --> gga and pp (to be override) are set automatically
	#  vdw-DFs: vdw-df, optpbe-vdw, optb88-vdw, optb86b-vdw, vdw-df2, beef-vdw
	#    --> luse_vdw and others are set automatically
	xc          = "beef-vdw"
	ivdw        = 0
	prec        = "normal"
	encut       = 400.0  # 213.0 or 400.0 or 500.0
	potim       = 0.10
	ibrion      = 2
	nsw         = 200
	nsw_neb     = 20
	nsw_dimer   = 1000
	nelmin      = 5
	nelm        = 50  # default:40
	ediff       = 1.0e-6
	ediffg      = -0.05
	kpts        = [3, 3, 1]
	ismear      = 1
	sigma       = 0.10
	vacuum      = 10.0  # for gas-phase molecules. surface vacuum is set by surf.py
	setups      = None
	ialgo       = 48  # normal=38, fast=58, veryfast=48
	lwave       = False
	lcharg      = False
	ispin       = 1
	#setups = {"O" : "_h"}

	npar = 18  # ito: 18 hokudai: 10
	nsim = npar

	# set lmaxmix
	lmaxmix = 2
	first_TM = ["Fe", "Ni"]
	for elem in first_TM:
		if elem in surf.get_chemical_symbols():
			lmaxmix = 4

	if xc == 'pw91':
		pp = 'potpaw_GGA'
	else:  # mostly PBE is OK
		pp = 'potpaw_PBE.54'

	# dipoole correction
	dipole = True
	if dipole:
		ldipol, idipol = True, 3

	method = xc
	basis  = ""
	label  = method

	# switch of ivdw for vdw-including xc
	if xc in ["vdw-df", "optpbe-vdw", "optb88-vdw", "optb86b-vdw", "vdw-df2", "beef-vdw"]:
		ivdw = 0

	# DFT+U
	DFTU = False
	if DFTU:
		ldau = "true"
		ldautype = 2
		ldau_luj = {'La': {'L': 3, 'U': 5.0, 'J': 0.0}, 'O': {'L': -1, 'U': 0.0, 'J': 0.0}}
		#ialgo = 38

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

	# ---
	calc_gas  = Vasp(label=dir, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw,
					 npar=npar, nsim=nsim, encut=encut, ismear=0, istart=0, setups=setups, sigma=sigma,
					 ialgo=ialgo, lwave=lwave, lcharg=lcharg, ibrion=ibrion, potim=potim, nsw=nsw,
					 ediff=ediff, ediffg=ediffg, kpts=[1, 1, 1], pp=pp, ldipol=ldipol, idipol=idipol)
	calc_surf = Vasp(label=dir, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw,
					 npar=npar, nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
					 ialgo=ialgo, lwave=lwave, lcharg=lcharg, ibrion=ibrion, potim=potim, nsw=nsw,
					 ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol, idipol=idipol)
	if SP:
		calc_gas_sp  = Vasp(label=dir, prec=prec, xc=xc_sp, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw_sp,
							npar=npar, nsim=nsim, encut=encut, ismear=0, istart=0, setups=setups, sigma=sigma,
							ialgo=ialgo_sp, lwave=lwave, lcharg=lcharg, ibrion=ibrion, potim=potim, nsw=nsw,
							ediff=ediff, ediffg=ediffg, kpts=[1, 1, 1], pp=pp, ldipol=ldipol, idipol=idipol)
		calc_surf_sp = Vasp(label=dir, prec=prec, xc=xc_sp, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw_sp,
					    	npar=npar, nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups,
							sigma=sigma, ialgo=ialgo_sp, lwave=lwave, lcharg=lcharg, ibrion=ibrion, potim=potim,
							nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, pp=pp, ldipol=ldipol, idipol=idipol)

	if any(IR):
		calc_ir = Vasp(label=irdir, prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw,
					   npar=npar, nsim=nsim, encut=encut, ismear=ismear, setups=setups,
					   sigma=sigma, ialgo=ialgo, lwave=lwave, lcharg=lcharg, ediff=ediff, kpts=kpts,
					   isim=0, lmono=True, dipol=dipol, idipol=idipol)
## --- EMT --- -> nothing to set
else:
	method = ""
	basis  = ""
	label  = "emt"
	fmax   = 0.10

if True in ZPE:
	label = label + "_ZPE"
if SP:
	label = label + "_SP"

barrierfile = reactionfile.split(".")[0] + "_Ea_" + label + ".txt"
deltaEfile  = "deltaE.txt"
fbarrier = open(barrierfile, 'w')
fdeltaE  = open(deltaEfile, 'w')
fbarrier.close()
fdeltaE.close()

print("calculator: %s, method: %s, basis: %s" % (calculator, method, basis))
#
# start calculation
#
for irxn in range(rxnst, rxned):
	fbarrier = open(barrierfile, 'a')
	fdeltaE  = open(deltaEfile,  'a')
	print("--- calculating elementary reaction No. %d ---" % irxn)

	energies_bothside = {"reactant": None, "product": None}
	#reac_en = np.array(range(len(r_ads[irxn])), dtype="f")
	#prod_en = np.array(range(len(p_ads[irxn])), dtype="f")
	dirs = {"reactant": None, "product": None}  # working directories. set by label in Vasp calculator.

	for side in ["reactant", "product"]:
		if side == "reactant":
			mol_set = r_ads.copy()
			sites   = r_site.copy()
			coefs   = r_coef.copy()
		elif side == "product":
			mol_set = p_ads.copy()
			sites   = p_site.copy()
			coefs   = p_coef.copy()

		energies = np.array(range(len(mol_set[irxn])), dtype="f")

		for imols, mols in enumerate(mol_set[irxn]):
			surf_tmp = surf.copy()

			for imol, mol in enumerate(mols):
				ads_height = ads_height0
				print("----- %s: molecule No. %d is %s -----" % (side, imol, mol))
				config = "normal"

				if 'surf' in mol:
					# surface itself
					mol, neutral, charge = read_charge(mol)
					tmp = surf
				elif 'def' in mol:
					# defected surface
					mol, neutral, charge = read_charge(mol)
					tmp = 'def'
				else:
					# some adsorbate is present
					mol, neutral, charge = read_charge(mol)

					# flip or rotate
					if "-SIDEy" in mol:
						mol = mol.replace("-SIDEy", "")
						tmp = methane[mol]
						tmp.rotate(90, 'y')
						config = "sidey"
					elif "-SIDEx" in mol:
						mol = mol.replace("-SIDEx", "")
						tmp = methane[mol]
						tmp.rotate(90, 'x')
						config = "sidex"
					elif "-FLIP" in mol:
						mol = mol.replace("-FLIP", "")
						tmp = methane[mol]
						tmp.rotate(180, 'y')
						config = "flip"
					elif "-TILT" in mol:
						mol = mol.replace("-TILT", "")
						tmp = methane[mol]
						tmp.rotate(45, 'y')
						config = "tilt"
					elif "-HIGH" in mol:
						mol = mol.replace("-HIGH", "")
						tmp = methane[mol]
						ads_height += 1.5
						config = "high"
					elif "-AIR" in mol:
						mol = mol.replace("-AIR", "")
						tmp = methane[mol]
						ads_height += 4.0
						config = "air"
					else:
						tmp = methane[mol]

				site = sites[irxn][imols][imol]

				try:
					site, site_pos = site.split(".")
				except:
					site_pos = 'x1y1'

				if site == 'gas':
					if 'surf' in mols:
						mol_type = 'surf'
					else:
						mol_type = 'gaseous'
				else:
					mol_type = 'adsorbed'
			
				if mol_type == 'adsorbed':
					offset = site_info[lattice][facet][site][site_pos]
					if len(offset) == 3:
						shift = offset[2] * surf.get_cell()[2][2]
						ads_height += shift
						offset = offset[0:2]

					offset = np.array(offset)
					#
					# wrap atoms to prevent adsorbate being on different cell
					#
					surf_tmp.translate([0, 0, 2])
					surf_tmp.wrap(pbc=[0, 0, 1])
					surf_tmp.translate([0, 0, -2])
					print("lattice:{0}, facet:{1}, site:{2}, site_pos:{3}, config:{4}"
										.format(lattice, facet, site, site_pos, config))
					#
					if tmp == 'def':
						defect = find_closest_atom(surf_tmp, offset=offset*offset_fac)
						del surf_tmp[len(surf_tmp.get_atomic_numbers())-1]
						del surf_tmp[defect]  # defect
						tmp = surf_tmp
					else:
						#
						# shift adsorbate molecule
						#
						xcom  = tmp.get_center_of_mass()[0]
						ycom  = tmp.get_center_of_mass()[1]
						shift = (xcom - tmp.positions[0, 0], ycom - tmp.positions[0, 1])
						z_shift = tmp.positions[:, 2].min()

						if tmp.get_chemical_formula()  == 'H':  # special attention to H
							ads_height = 1.0

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
				print("found in database")
				past = tmpdb.get(name=formula + site + site_pos + config)
				first_time = False
			except:
				print("first time")
				first_time = True
			#else:
			#	if site == past.data.site:
			#		if site_pos == past.data.site_pos:
			#			if config == past.data.config:
			#				if len(mols) == 1:
			#					print("already calculated")
			#					tmp = tmpdb.get_atoms(id=past.id)
			#					first_time = False

			magmom = tmp.get_initial_magnetic_moments()
			natom  = len(tmp.get_atomic_numbers())
			coef   = coefs[irxn][imols]

			# spin switch
			if int(math.ceil(sum(magmom))) != 0:
				print("has unpaired electron")
				ispin = 2
			else:
				ispin = 1

			if natom == 1:
				ZPE[side] = False
				IR[side]  = False
			#
			# set label
			#
			dir = label + "_rxn" + str(irxn).zfill(3) + "_"
			if mol_type == "surf":
				dir += surf_name
			elif mol_type == "adsorbed":
				dir += "-".join(mols) + site + "_" + surf_name
			else:
				dir += "-".join(mols)

			if not os.path.isdir(traj_dir):
				os.makedirs(traj_dir)
			traj = traj_dir + "/" + dir + "_" + side[0:4] + ".traj"
			#
			# branch compurational setting by gas or not
			#
			if site == "gas" and "surf" not in mols:
				#
 				# gas-phase molecule
				#
				gas_mol = True
				if calculator == "vasp":
					cell = np.array([1, 1, 1])
					cell = vacuum*cell
					tmp.set_cell(cell)
					tmp.center()
			else:  # surface
				gas_mol = False
			#
			# set calculator
			#
			if calculator == "gaussian":
				tmp.calc = calc_gas
				opt = BFGS(tmp, trajectory=traj)
				opt.run(fmax=fmax, steps=maxoptsteps)
				if SP:
					spdir = dir + "_sp"
					tmp.calc = Gaussian(label=spdir, method=method_sp, basis=basis,
										force=None, nprocs=nprocs, mem=mem)
			elif calculator == "vasp":
				if gas_mol:
					#
					# gas-phase molecules
					#
					tmp.calc = calc_gas
					if not neutral:
						nelect = get_number_of_valence_electrons(tmp)
						nelect = nelect - charge
						tmp.calc.int_params["nelect"] = nelect
						tmp.calc.bool_params["lmono"] = "True"
				else:
					# surface
					tmp.calc = calc_surf
					if not neutral:
						nelect = get_number_of_valence_electrons(tmp)
						nelect = nelect - charge
						tmp.calc.int_params["nelect"] = nelect
						tmp.calc.bool_params["lmono"] = "True"
					if DFTU:
						tmp.calc.bool_params["ldau"] = ldau
						tmp.calc.int_params["ldautype"] = ldautype
						tmp.calc.dict_params["ldau_luj"] = ldau_luj

			elif "emt" in calculator:
				tmp.calc = EMT()
				opt = BFGS(tmp, trajectory=traj)
				opt.run(fmax=fmax, steps=maxoptsteps)

			en = tmp.get_potential_energy()

			# single point
			if SP:
				if mol_type == 'gaseous':
					tmp.calc = calc_gas_sp
				else:
					tmp.calc = calc_surf_sp
				en = tmp.get_potential_energy()

			if calculator == "vasp":
				xmlfile = "vasprun_" + dir + ".xml"
				contcar = "CONTCAR_" + dir
				reac_contcar = contcar
				os.system('cp vasprun.xml %s >& /dev/null' % xmlfile)
				os.system('cp CONTCAR %s >& /dev/null' % contcar)
				os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHGCAR CHG WAVECAR REPORT >& /dev/null')

			if mol_type != 'surf' and (ZPE[side] or IR[side]):  # surf--nothing to do with vibration
				# fix atoms for vibrations
				if mol_type == 'adsorbed':
					c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1 or atom.tag == 2])
					tmp.set_constraint(c)
				if ZPE[side]:
					# do ZPE calculation
					zpedir = dir + "_ZPE"
					if not os.path.isdir(zpedir):
						os.makedirs(zpedir)
					vib = Vibrations(tmp, delta=0.01, name=zpedir + "/" + "ZPE")  # delta = 0.01 is default
					vib.run()
					hnu = vib.get_energies()
					zpe = vib.get_zero_point_energy()
					en = en + zpe
					os.system("rm -r zpedir")

				if IR[side]:
					# setting for IR calculation
					irdir = dir + "_IR"
					if not os.path.isdir(irdir):
						os.makedirs(irdir)

					tmp.calc  = calc_ir
					tmp.calc.list_float_params[dipol] = tmp.get_center_of_mass(scaled=True)

					vib = Infrared(tmp, delta=0.01, name=irdir+"/"+"IR")
					vib.run()
					vib.write_spectra(out=irdir + "_spectra.dat", start=1000, end=4000)
					vib.write_mode()
					vib.write_jmol()

					# include ZPE
					zpe = vib.get_zero_point_energy()
					en = en + zpe
					os.system("rm *.pckl")

				energies[imols] = en
			else:
				energies[imols] = en

			energies[imols] = coef * energies[imols]

			dirs[side] = dir  # needed when doing TS calc

			# recording to database
			if first_time:
				id = tmpdb.reserve(name=formula + site + site_pos + config)
				if id is None:  # somebody is writing to db
					continue
				else:
					tmpdb.write(tmp, name=formula + site + site_pos + config, id=id,
								data={'site': site, 'site_pos': site_pos, 'config': config})

		energies_bothside[side] = energies

	deltaE = np.sum(energies_bothside["product"]) - np.sum(energies_bothside["reactant"])

	print("deltaE = %5.3f" % deltaE)

	#
	# end loop over side
	#

	#
	# TS calc
	#
	if TS:
		# make different directory and go there
		os.system('rm -rf tsdir')  # remove old one

		if not os.path.exists("tsdir"):
			os.makedirs("tsdir")
		os.chdir("tsdir")

		# copy reactant and product POSCAR as POSCAR_reac and POSCAR_prod
		contcar1 = "../" + dirs["reactant"] + "/CONTCAR"
		contcar2 = "../" + dirs["product"]  + "/CONTCAR"
		os.system('cp %s ./POSCAR_reac' % contcar1)
		os.system('cp %s ./POSCAR_prod' % contcar2)

		# Swith atomic index in order to make POSCAR_reac and POSCAR_prod close.
		# Different ordering cause bad NEB images.
		atom1 = read('POSCAR_reac')
		atom2 = read('POSCAR_prod')
		newatom1 = make_it_closer_by_exchange(atom1, atom2, thre=100.0)  # atom1 is exchanged
		write('POSCAR_reac', newatom1)  # do not sort because POTCAR does not follow
		write('POSCAR_prod', atom2)
		#
		# do "nebmake.pl"
		#
		nebmake = vtst + "nebmake.pl"

		os.system('%s POSCAR_reac POSCAR_prod %d >& /dev/null' % (nebmake, nimages))

		outcar1  = "../" + dirs["reactant"] + "/OUTCAR"
		outcar2  = "../" + dirs["product"] + "/OUTCAR"
		# copy reactant and product CONTCAR and OUTCAR to 00 and 0#IMAGE's POSCAR and OUTCAR
		os.system('cp %s 00'          % outcar1)
		os.system('cp %s 00/POSCAR'   % contcar1)
		os.system('cp %s %02d'        % (outcar2,  nimages+1))
		os.system('cp %s %02d/POSCAR' % (contcar2, nimages+1))

		# normal NEB
		if neutral:
			tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
							nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
							ialgo=ialgo, lwave=lwave, lcharg=lcharg,
		 					ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts,
		 					images=nimages, spring=-5.0, lclimb=False, iopt=7, maxmove=0.10, pp=pp,
							ldipol=ldipol, idipol=idipol)
		else:
			nelect = get_number_of_valence_electrons(tmp)
			nelect = nelect - charge
			tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
							nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
							ialgo=ialgo, lwave=lwave, lcharg=lcharg,
		 					ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts,
		 					images=nimages, spring=-5.0, lclimb=False, iopt=7, maxmove=0.10, nelect=nelect,
							lmono="true", pp=pp, ldipol=ldipol, idipol=idipol)

		print("----------- doing NEB calculation with images=%d ----------" % nimages)
		tmp.get_potential_energy()
		print("----------- normal NEB done -----------")
		neb_copy_contcar_to_poscar(nimages)

		# CI-NEB
		if CI:
			if neutral:
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
								nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
								ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts,
								images=nimages, spring=-5.0, lclimb=True, iopt=7, maxmove=0.10, pp=pp,
								ldipol=ldipol, idipol=idipol)
			else:
				nelect = get_number_of_valence_electrons(tmp)
				nelect = nelect - charge
				tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
								nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
								ialgo=ialgo, lwave=lwave, lcharg=lcharg,
								ibrion=3, potim=0, nsw=nsw_neb, ediff=ediff, ediffg=ediffg, kpts=kpts,
								images=nimages, spring=-5.0, lclimb=True, iopt=7, maxmove=0.10, nelect=nelect,
								lmono="true", pp=pp, ldipol=ldipol, idipol=idipol)
			print("---------- doing CI-NEB calculation with images=%d -------" % nimages)
			tmp.get_potential_energy()
			print("----------- CI NEB done -----------")
			neb_copy_contcar_to_poscar(nimages)

		nebresults = vtst + "nebresults.pl"
		neb2dim    = vtst + "neb2dim.pl"
		os.system('%s >& /dev/null' % nebresults)
		os.system('%s >& /dev/null' % neb2dim)
		os.chdir("dim")

		# print("dimer made"); quit()

		# dimer method
		if neutral:
			tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
							nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
							ialgo=ialgo, lwave=lwave, lcharg=lcharg,
							nsw=nsw_dimer, ediff=ediff*0.1, ediffg=ediffg*0.5, kpts=kpts,
							iopt=2, maxmove=0.20, dfnmax=1.0, ichain=2, pp=pp, ldipol=ldipol, idipol=idipol)
		else:
			nelect = get_number_of_valence_electrons(tmp)
			nelect = nelect - charge
			tmp.calc = Vasp(prec=prec, xc=xc, ispin=ispin, nelm=nelm, nelmin=nelmin, ivdw=ivdw, npar=npar,
							nsim=nsim, encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma,
							ialgo=ialgo, lwave=lwave, lcharg=lcharg,
							nsw=nsw_dimer, ediff=ediff*0.1, ediffg=ediffg*0.5, kpts=kpts,
							iopt=2, maxmove=0.20, dfnmax=1.0, ichain=2, nelect=nelect, lmono="true", pp=pp,
							ldipol=ldipol, idipol=idipol)

		print("----------- doing dimer method TS optimization -----------")
		TS_en = tmp.get_potential_energy()
		print("----------- dimer done -----------")

		Ea = TS_en - reac_en
		print("Ea = %5.3f" % Ea)

		# TS end

	#
	# writing reaction
	#
	string = ""
	for side in ["reactant", "product"]:
		if side == "reactant":
			mols = r_ads
			sites = r_site
		elif side == "product":
			mols = p_ads
			sites = p_site

		for imol, mol in enumerate(mols[irxn]):
			mol_join = '-'.join(mol)
			string = string + "{0}_{1}".format(mol_join, '-'.join(sites[irxn][imol]))
			if imol != len(mols[irxn])-1:
				string = string + " + "

		if side == "reactant":
			string += " --> "

	fbarrier.write('{0:>3d} {1:<70s}'.format(irxn, string))

	Ea_for  = deltaE
	Ea_rev  = -deltaE
	Ea = [Ea_for, Ea_rev]

	fbarrier.write('{0:>14.8f} {1:>14.8f}\n'.format(Ea_for, Ea_rev))
	fbarrier.close()

	fdeltaE.write('{0:>3d} {1:>14.8f} {2:>14.8f}\n'.format(irxn, Ea_for, Ea_rev))
	fdeltaE.close()
	#
	# end loop over reaction
	#

remove_parentheses(barrierfile)
