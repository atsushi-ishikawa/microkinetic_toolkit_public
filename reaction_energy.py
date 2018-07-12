import numpy as np
import os, sys, re, json
from reaction_tools import *
from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.collections import methane
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.db import connect
from ase.io import read
from ase.build import add_adsorbate
from ase.visualize import view
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
# if surface present, provide surface file
# in ase.db form
#
surface = True

if surface:
	db = connect('surf.db')
	surf      = db.get_atoms(id=1)
	lattice   = db.get(id=1).data.lattice
	facet     = db.get(id=1).data.facet
	surf_name = db.get(id=1).data.formula

	# load site information
	f = open('site_info.json','r')
	site_info = json.load(f)
	f.close()

# fix atoms
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1])
surf.set_constraint(c)

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

rxn_num = get_number_of_reaction(reactionfile)
Ea = np.array(2, dtype="f")

## --- parameters
ZPE = False
SP  = False
maxoptsteps = 200
ads_height0 = 1.4
ads_pos = (0.0, 0.0)
# whether to do single point after optimization
# at different computational level

## --- Gaussian ---
if "gau" in calculator:
	method = "b3lyp"
	basis  = "6-31G(d)" # do not use aesterisk for polarization func
	if SP:
		method_sp = "ccsd(t)"
	basis_name = re.sub("\(", "", basis)
	basis_name = re.sub("\)", "", basis_name)
	basis_name = re.sub(",",  "", basis_name)
	label = method + "_" + basis_name

## --- VASP ---
elif "vasp" in calculator:
	xc          = "rpbe"
	prec        = "low"
	encut       = 350.0 # 213.0 or 400.0 or 500.0
	potim       = 0.10
	nsw         = 10
	nelmin      = 5
	ediff       = 1.0e-4
	ediffg      = -0.1
	kpts_surf   = [1, 1, 1]
	ismear_surf = 1
	sigma_surf  = 0.20
	vacuum      = 10.0 # for gas-phase molecules. surface vacuum is set by surf.py
	setups      = None
	ivdw        = 12
	ialgo       = 48 # normal=38, veryfast=48
	npar        = 12
	nsim        = npar
	#setups = {"O" : "_h"}

	method = xc
	basis  = ""
	label  = method

	# DFT+U
	DFTU = True
	if DFTU:
		ldau = "true"
		ldautype = 2
		ldau_luj = { 'La':{'L':3, 'U':3.5, 'J':0.0}, 'O':{'L':-1, 'U':0.0, 'J':0.0} }
		ialgo = 38
	# charge
	neutral = True

## --- EMT --- -> nothing to set
elif "emt" in calculator:
	method = ""
	basis  = ""
	label  = ""

if ZPE:
	label = label + "ZPE"
if SP:
	label = label + "SP"

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
	reac_A  = np.array(range(len(r_ads[irxn])),dtype="f")
	prod_A  = np.array(range(len(r_ads[irxn])),dtype="f")

	#
	# reactants
	#
	for imols, mols in enumerate(r_ads[irxn]):
		surf_tmp = surf.copy()
		for imol, mol in enumerate(mols):
			print "----- reactant: molecule No.", imol, " is ", mol, "-----"

			if 'surf' in mol:
				mol,neutral,charge = read_charge(mol)
				tmp = surf
			elif 'def' in mol:
				mol,neutral,charge = read_charge(mol)
				tmp = 'def'
			else:
				mol,neutral,charge = read_charge(mol)

				# flip or rotate
				if "-SIDE" in mol:
					mol = mol.replace("-SIDE","")
					tmp = methane[mol]
					tmp.rotate(90,'y')
				elif "-FLIP" in mol:
					mol = mol.replace("-FLIP","")
					tmp = methane[mol]
					tmp.rotate(180,'y')
				else:
					tmp = methane[mol]

			site = r_site[irxn][imols][imol]

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			if site != 'gas':
				offset = site_info[lattice][facet][site][site_pos]
				offset = np.array(offset)*(3.0/4.0) # MgO only
				# wrap atoms to prevent adsorbate being on different cell
				surf_tmp.translate([0,0,2])
				surf_tmp.wrap(pbc=[0,0,1])
				surf_tmp.translate([0,0,-2])
				print("lattice:{0}, facet:{1}, site:{2}, site_pos:{3}\n".format(lattice,facet,site,site_pos))
				#
				if tmp == 'def':
					defect = find_closest_atom(surf_tmp,offset=offset)
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
						ads_height = 0.9

					ads_height = ads_height0 - z_shift
					ads_pos = (-shift[0], -shift[1])
					add_adsorbate(surf_tmp, tmp, ads_height, position=ads_pos, offset=offset)
					tmp = surf_tmp
		del surf_tmp
		#
		# end adsorbing molecule
		#
		magmom = tmp.get_initial_magnetic_moments()
		natom  = len(tmp.get_atomic_numbers())
		coef   = r_coef[irxn][imols]
		#
		# set label
		#
		r_label = label + "_rxn" + str(irxn) + "_" + "-".join(mols) + "_" + site
		if site != 'gas':
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
			ismear  = ismear_surf # Methfessel-Paxton
			sigma   = sigma_surf
			kpts    = kpts_surf
			gas_mol = False
		#
		# set calculator
		#
		if "gau" in calculator:
			tmp.calc = Gaussian(label=r_label, method=method, basis=basis, scf="xqc")
			opt = BFGS(tmp, trajectory=r_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)
			if SP:
				r_label = r_label + "_sp"
				tmp.calc = Gaussian(label=r_label, method=method_sp, basis=basis, force=None)
		elif "vasp" in calculator:
			if gas_mol:
				#
				# gas-phase molecules
				#
				if neutral:
		 			tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
									ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts )
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
		 			tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
									ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
									nelect=nelect, lmono="true" )
			else:
				#
				# surface
				#
				if neutral:
					if DFTU:
		 				tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj ) # DFTU
					else:
			 			tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts) # normal
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
					if DFTU:
		 				tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true", ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj ) # charge + DFTU
					else:
		 				tmp.calc = Vasp(output_template=r_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true") # charge

		elif "emt" in calculator:
			tmp.calc = EMT()
			opt = BFGS(tmp, trajectory=r_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)

		en = tmp.get_potential_energy()

		if "vasp" in calculator:
			xmlfile = "vasprun_" + r_label + ".xml"
			os.system('cp vasprun.xml %s' % xmlfile)
			os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHGCAR CHG WAVECAR REPORT')
		if ZPE == True and natom != 1:
			vib = Vibrations(tmp)
			vib.run()
			hnu = vib.get_energies()
			zpe = vib.get_zero_point_energy()
			reac_en[imols] = en + zpe
			os.system("rm vib.*")
		else:
			reac_en[imols] = en

		reac_en[imols] = coef * reac_en[imols]
	#
	# products
	#
	for imols, mols in enumerate(p_ads[irxn]):
		surf_tmp = surf.copy()
		for imol, mol in enumerate(mols):
			print "----- product: molecule No.", imol, " is ", mol, "-----"

			if 'surf' in mol:
				mol,neutral,charge = read_charge(mol)
				tmp = surf
			elif 'def' in mol:
				mol,neutral,charge = read_charge(mol)
				tmp = 'def'
			else:
				mol,neutral,charge = read_charge(mol)

				# flip or rotate
				if "-SIDE" in mol:
					mol = mol.replace("-SIDE","")
					tmp = methane[mol]
					tmp.rotate(90,'y')
				elif "-FLIP" in mol:
					mol = mol.replace("-FLIP","")
					tmp = methane[mol]
					tmp.rotate(180,'y')
				else:
					tmp = methane[mol]

			site = p_site[irxn][imols][imol]

			try:
				site,site_pos = site.split(".")
			except:
				site_pos = 'x1y1'

			if site != 'gas':
				offset = site_info[lattice][facet][site][site_pos]
				offset = np.array(offset)*(3.0/4.0) # MgO only
				# wrap atoms to prevent adsorbate being on different cell
				surf_tmp.translate([0,0,2])
				surf_tmp.wrap(pbc=[0,0,1])
				surf_tmp.translate([0,0,-2])
				print("lattice:{0}, facet:{1}, site:{2}, site_pos:{3}\n".format(lattice,facet,site,site_pos))
				#
				if tmp == 'def':
					defect = find_closest_atom(surf_tmp,offset=offset)
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
						ads_height = 0.9

					ads_height = ads_height0 - z_shift
					ads_pos = (-shift[0], -shift[1])
					add_adsorbate(surf_tmp, tmp, ads_height, position=ads_pos, offset=offset)
					tmp = surf_tmp
		del surf_tmp
		#
		# end adsorbing molecule
		#
		magmom = tmp.get_initial_magnetic_moments()
		natom  = len(tmp.get_atomic_numbers())
		coef   = p_coef[irxn][imols]
		#
		# set label
		#
		p_label = label + "_rxn" + str(irxn) + "_" + "-".join(mols) + "_" + site
		if site != 'gas':
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
			tmp.calc = Gaussian(label=p_label, method=method, basis=basis, scf="xqc")
			opt = BFGS(tmp, trajectory=p_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)
			if SP:
				p_label = p_label + "_sp"
				tmp.calc = Gaussian(label=p_label, method=method_sp, basis=basis, force=None)
		elif "vasp" in calculator:
			if gas_mol:
				#
				# gas-phase molecules
				#
				if neutral:
		 			tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
									ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts )
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
		 			tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
									encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
									ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts,
									nelect=nelect, lmono="true" )
			else:
				#
				# surface
				#
				if neutral:
					if DFTU:
		 				tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj ) # DFTU
					else:
			 			tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts) # normal
				else:
					nelect = get_number_of_valence_electrons(tmp)
					nelect = nelect - charge
					if DFTU:
		 				tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true", ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj ) # charge + DFTU
					else:
		 				tmp.calc = Vasp(output_template=p_label, prec=prec, xc=xc, ispin=2, nelmin=nelmin, ivdw=ivdw, npar=npar, nsim=nsim,
										encut=encut, ismear=ismear, istart=0, setups=setups, sigma=sigma, ialgo=ialgo,
										ibrion=2, potim=potim, nsw=nsw, ediff=ediff, ediffg=ediffg, kpts=kpts, 
										nelect=nelect, lmono="true") # charge

		elif "emt" in calculator:
			tmp.calc = EMT()
			opt = BFGS(tmp, trajectory=p_traj)
			opt.run(fmax=0.05, steps=maxoptsteps)

		en = tmp.get_potential_energy()

		if "vasp" in calculator:
			xmlfile = 'vasprun_' + p_label + '.xml'
			os.system('cp vasprun.xml %s' % xmlfile)
			os.system('rm PCDAT XDATCAR EIGENVAL OSZICAR IBZKPT CHGCAR CHG WAVECAR REPORT')
		if ZPE == True and natom != 1:
			vib = Vibrations(tmp)
			vib.run()
			hnu = vib.get_energies()
			zpe = vib.get_zero_point_energy()
			prod_en[imols] = en + zpe
			os.system('rm vib.*')
		else:
			prod_en[imols] = en

		prod_en[imols] = coef * prod_en[imols]

	deltaE = np.sum(prod_en) - np.sum(reac_en)
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

