from microkinetic_toolkit.reaction import Reaction
import os
import numpy as np
import pandas as pd

class Reactions:
	"""
	Set of elementary reactions.
	"""
	def __init__(self, reaction_list):
		self.reaction_list = reaction_list

	def __getitem__(self, index):
		return self.reaction_list[index]

	def __len__(self):
		return len(self.reaction_list)

	def to_tdb(self, db_file, update=False):
		tdb = TinyDB(db_file)
		for reaction in self.reaction_list:
			if update:
				reaction.update_tdb(tdb)
			else:
				reaction.to_tdb(tdb)

	def to_openfoam(self, file):
		"""
		Generate openFOAM input file.
		Args:
			file: name of the generated openFOAM input file
		"""
		with open(file, "w") as write:
			species_info = self._get_openfoam_species_info()
			write.write(species_info)
			for info in self._get_openfoam_reactions_info():
				write.write(info)

	def _get_openfoam_species_info(self):
		ini = "species\n(\n"
		fin = "\n)\n\n"
		species = "\n".join([str(sp) for sp in self.get_unique_species()])
		species = textwrap.indent(species, prefix="    ")
		species_info = ini + species + fin
		return species_info

	def get_unique_species(self):
		"""
		Get unique chemical species.
		Returns:
			list of string
		"""
		species_set = set([])
		for reaction in self.reaction_list:
			species_set.update(reaction.unique_species)
		return list(species_set)

	def _get_openfoam_reactions_info(self):
		ini = "reaction\n{\n"
		yield ini
		for react_info in self._get_each_react_info():
			react_info = textwrap.indent(react_info, prefix="   ")
			yield react_info
		fin = "}\n"
		yield fin

	def _get_each_react_info(self):
		for reaction in self.reaction_list:
			paramdict = reaction.to_openfoam_paramdict()
			ini = "TestReaction{}\n".format(reaction._reaction_id)
			ini += "{\n"
			conts = "type     {0[type]};\n" \
					"reaction {0[reaction]};\n" \
					"A        {0[A]};\n" \
					"beta     {0[beta]};\n" \
					"Ta       {0[Ta]};\n"
			conts = conts.format(paramdict)
			conts = textwrap.indent(conts, prefix="    ")
			fin = "}\n"

			react_info = ini + conts + fin
			yield react_info

	def _generate_reactions_dict(self):
		for reaction in self.reaction_list:
			ddict = reaction.to_dict()
			yield ddict

	def to_csv(self, file):
		"""
		Generate csv file containing elemtary reactions.
		Args:
			file: csv file name
		"""
		df = DataFrame(self._generate_reactions_dict())
		df.to_csv(file)

	@classmethod
	def from_csv(cls, csv_file):
		"""
		Read elementary reactions from CSV.
		Args:
			csv_file: CSV file with elementary reactions
		Returns:
			Reactions
		"""
		df = pd.read_csv(csv_file, index_col=0)
		reaction_list = []
		for i, row in df.iterrows():
			ddict = row.to_dict()
			reaction = Reaction.from_dict(ddict)
			reaction_list.append(reaction)
		return cls(reaction_list)

	def adsorbate_on_surface(self, surface=None):
		"""
		Make surface with adsorbates.
		Args:
			surface:
		Returns:
		"""
		for reaction in self.reaction_list:
			reac, prod = reaction.adsorbate_on_surface(surface)
		pass

	def get_reaction_energies(self, surface=None, method=None):
		"""
		Calculate the reaction energies (deltaE) for all the elementary reactions.
		Args:
			surface: Atoms
			method: emt
		Returns:
			deltaEs: numpy array
		"""
		deltaEs = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			deltaEs[i] = reaction.get_reaction_energy(surface=surface, method=method)
		return deltaEs

	def get_entropy_differences(self):
		"""
		Calculate the entropy difference (deltaS, in eV/K) for all the elementary reactions.
		Returns:
			deltaSs: numpy array
		"""
		deltaSs = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			deltaSs[i] = reaction.get_entropy_difference()
		return deltaSs

	def get_rate_constants(self, deltaEs=None, T=300.0):
		"""
		Calculate rate constants for all the elementary reactions.
		Args:
			deltaEs: reaction energies [eV]
			T: temperature [K]
		Returns:
			ks: rate constants (numpy array)
		"""
		ks = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			index  = reaction._reaction_id
			deltaE = deltaEs[index]
			ks[i]  = reaction.get_rate_constant(deltaE, T)
		return ks

	def calculate_volume_and_entropy(self):
		"""
		Calculate volume and entropy for all the molecules in the current reaction set.
		Returns:
		"""
		from ase import Atoms, Atom
		from ase.calculators.gaussian import Gaussian
		from ase.db import connect
		from ase.optimize import BFGS
		from ase.vibrations import Infrared
		#
		# calculate reaction energy
		# molecule's data should be stored in "methane.json"
		#
		calculator = "gau"
		calculator = calculator.lower()

		oldjson = "tmp.json"
		newjson = "tmp_new.json"

		db1 = connect(oldjson)
		db2 = connect(newjson)
		#
		add_volume = True
		add_entropy = True
		#
		###
		## --- Gaussian ---
		if "gau" in calculator:
			method = "b3lyp"
			basis = "6-31G*"
		else:
			raise RuntimeError('molecular volume only implemented in Gaussian')
		#
		# add volume and entropy information for all the molecules in database
		#
		Nmol = len(db1)
		for imol in range(Nmol):  # add 10 for safety
			try:
				id = imol + 1
				mol = db1.get(id=id).name
				tmp = db1.get_atoms(id=id)
				print("adding volume and entropy", mol)

				# do geometry optimization first
				tmp.calc = Gaussian(method=method, basis=basis)
				BFGS(tmp).run(fmax=0.05)

				# volume and molecular total entropy calculation
				if add_volume and not add_entropy:
					tmp.calc = Gaussian(method=method, basis=basis, volume='tight')
					tmp.get_potential_energy()
					vol = tmp.get_molecular_volume()
				if not add_volume and add_entropy:
					tmp.calc = Gaussian(method=method, basis=basis, force=None, freq='noraman')
					tmp.get_potential_energy()
					entropy = tmp.get_molecular_entropy()
				if add_volume and add_entropy:
					tmp.calc = Gaussian(method=method, basis=basis, volume='tight', force=None, freq='noraman')
					tmp.get_potential_energy()
					vol = tmp.get_molecular_volume()
					entropy = tmp.get_molecular_entropy()

				# look for magmom
				try:
					magmom = tmp.get_magnetic_moments()
				except:
					magmom = tmp.get_initial_magnetic_moments()

				tmp.set_initial_magnetic_moments(magmom)

				# now write to database
				if add_volume and add_entropy:
					db2.write(tmp, key_value_pairs={'name': mol,
													'molecular_volume': vol, 'molecular_entropy': entropy})
				elif add_volume and not add_entropy:
					db2.write(tmp, key_value_pairs={'name': mol,
													'molecular_volume': vol})
				elif not add_volume and add_entropy:
					db2.write(tmp, key_value_pairs={'name': mol,
													'molecular_entropy': entropy})
			except:
				print("has nothing --- go to next")

		return None

	# microkinetics
	def do_microkinetics(self, rateconstants=None, T=300.0, P=1.0):
		"""
		Do microkinetic analysis.
		Args:
			rateconstants: rate constants in forward direction.
			T: temperature [K]
			P: total pressure [bar]
		Returns:
			None
		"""
		if rateconstants is None:
			print("rate constant not found")
			exit(1)

		odefile = "tmpode.py"
		self.make_rate_equation(odefile=odefile)
		self.solve_rate_equation(odefile=odefile, rateconstants=rateconstants, T=T, P=P)
		return None

	def make_rate_equation(self, odefile=None):
		"""
		Make rate equation file
		Args:
			odefile: filename to write ODE equations.
		Returns:
			None
		"""
		import microkinetic_toolkit.reaction
		if odefile is None:
			raise ValueError("ODE file not found")

		#(r_ads, r_site, r_coef, p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

		# r_ads and p_ads are species list of ALL the elementary reactions.
		# e.g. if inputfile contains
		#  (1) A1 + B1 --> C1 + D1
		#  (2) A2 + B2 --> C2 + D2
		# it gives
		# r_ads = [ [['A1'],['B1']] , [['A2'],['B2']] ]
		# p_ads = [ [['C1'],['D1']] , [['C2'],['D2']] ]

		dirname = os.path.dirname(microkinetic_toolkit.reaction.__file__)
		fout = open(dirname + "/" + odefile, "w")

		fout.write('import numpy as np')
		fout.write("\n\n")
		fout.write('def func(t, c, kfor, Kc, T, sden, area, Vr, ngas, ncomp):')
		fout.write("\n\n")
		#
		# template - start
		#
		lines = [
		"\tR = 8.314e-3  # kJ/mol/K\n",
		"\tkrev = kfor / Kc\n",
		"\ttheta = c[0:ncomp]\n",
		"\ttheta = theta * sden\n"
		]
		fout.writelines(lines)
		#
		# template - end
		#
		nspecies = len(self.get_unique_species())
		fout.write("\trate = np.zeros(" + str(nspecies) + ")\n\n")

		dict1 = {}
		dict2 = {}
		for irxn in range(len(self.reaction_list)):
			rxn_idx  = str(irxn)
			reaction = self[irxn]

			# hash-tag based reactant and product species list FOR THIS REACTION
			list_r = []
			list_p = []
			#
			# making dict1, a dictionary with hash of species-number and molecule
			#
			for side in ["reactant", "product"]:
				if side == "reactant":
					terms = reaction.reactants
					list = list_r
				elif side == "product":
					terms = reaction.products
					list = list_p
				else:
					print("error asdf")
					exit(1)

				for term in terms:
					spe, site = term[1], term[2]
					if site != 'gas':
						spe += "_surf"

					spe_num = self.get_unique_species().index(spe)
					list.append(spe_num)
					dict1[spe_num] = spe
			# done for dict1

			for side in ["reactant", "product"]:
				for direction in ["forward", "reverse"]:
					if side == "reactant" and direction == "forward":
						mol_list1 = reaction.reactants
						add_to = list_r
						mol_list2 = reaction.reactants  # list corresponding to add_to
						coefs = [i[0] for i in mol_list2]
						term = "kfor[" + rxn_idx + "]"
						sign = " - "
					elif side == "reactant" and direction == "reverse":
						mol_list1 = reaction.products
						add_to = list_r
						mol_list2 = reaction.reactants
						coefs = [i[0] for i in mol_list2]
						term = "krev[" + rxn_idx + "]"
						sign = " + "
					elif side == "product" and direction == "forward":
						mol_list1 = reaction.reactants
						add_to = list_p
						mol_list2 = reaction.products
						coefs = [i[0] for i in mol_list2]
						term = "kfor[" + rxn_idx + "]"
						sign = " + "
					elif side == "product" and direction == "reverse":
						mol_list1 = reaction.products
						add_to = list_p
						mol_list2 = reaction.products
						coefs = [i[0] for i in mol_list2]
						term = "krev[" + rxn_idx + "]"
						sign = " - "

					# making single term
					for mol in mol_list1:
						coef, spe, site = mol[0], mol[1], mol[2]

						if site != "gas":
							spe += "_surf"

						spe_num = self.get_unique_species().index(spe)
						if site == "gas":
							if spe == "surf":  # bare surface
								theta = "theta[" + str(spe_num) + "]"
							else:  # gas-phase molecule
								theta = "c[" + str(spe_num) + "]"
						else:  # adsorbed species
							theta = "theta[" + str(spe_num) + "]"

						power = coef
						if power != 1:
							theta = theta + "**" + str(power)

						term = term + "*" + theta

					for mem in add_to:
						if dict1[mem] == "surf":
							continue  # bare surface ... skip

						# CHECK
						coef = 0
						for imol, mol in enumerate(mol_list2):
							spe = mol[1]
							adsorbate = dict1[mem].split("_")[0]
							if spe == adsorbate:
								coef    = coefs[imol]

						if coef == 0:
							print("something wrong at coef 1")
							exit()

						sto_coef = str(float(coef))

						if mem in dict2:
							dict2[mem] = dict2[mem] + sign + sto_coef + "*" + term  # NEGATIVE
						else:
							dict2[mem] = sign + sto_coef + "*" + term

						if "theta" in dict2[mem]:
							dict2[mem] = dict2[mem] + "*" + "(area/Vr)"

		# vacancy site
		if 'surf' in dict1.values():  # only when surface is involved
			tmp = ""
			for imol, mol in enumerate(dict1):
				comp = dict1[imol]
				if 'surf' in comp and comp != 'surf':
					tmp = tmp + " -rate[" + str(imol) + "]"

			dict2[len(dict2)] = tmp

		comment = "\n\t# species --- "

		for imol, mol in enumerate(dict2):
			fout.write("\trate[{0}] ={1}  # {2}\n".format(imol, dict2[imol], dict1[imol]))
			comment += "%s = %s " % (imol, dict1[imol])
		comment += "\n"

		fout.write(comment)
		#
		# tempelate - start
		#
		lines = [
		"\tif ncomp > ngas:\n",
		"\t\trate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface\n"
		]
		#
		# tempelate - end
		#
		fout.writelines(lines)

		fout.write("\treturn rate\n")
		fout.write("\n")
		fout.close()

		return None

	def solve_rate_equation(self, rateconstants=None, odefile=None, T=300.0, P=1.0):
		"""
		Solve rate equations.
		Args:
			rateconstants: rate constant in forward direction
			odefile: ODE file
			T: temperature [K]
			P: total pressure [bar]
		Returns:
			None
		"""
		import numpy as np
		from scipy.integrate import solve_ivp
		import matplotlib.pyplot as plt
		import pickle
		from microkinetic_toolkit.tmpode import func

		if rateconstants is None:
			raise ValueError("rate constants not found")
		if odefile is None:
			raise ValueError("ODE file not found")

		np.set_printoptions(precision=3, linewidth=100)

		# constants
		R = 8.314 * 1.0e-3  # kJ/mol/K
		eVtokJ = 96.487

		### parameters
		# TODO: make it class properties
		Pin = P*1e5  # inlet pressure [Pascal]
		v0  = 1e-5   # vol. flowrate [m^3/sec]. 1 [m^2/sec] = 1.0e6 [mL/sec] = 6.0e7 [mL/min]

		sden  = 1.0e-05  # site density [mol/m^2]
		w_cat = 1.0e-3  # catalyst weight [kg]
		area  = 1000 * w_cat  # surface area. [m^2/kg] (e.g. BET) * [kg] --> [m^2]

		# Bronsted-Evans-Polanyi rule ... in eV unit
		alpha = 1.0
		beta  = 1.0

		phi = 0.5  # porosity
		rho_b = 1.0e3  # density of catalyst [kg/m^3]. typical is 1.0 g/cm^3 = 1.0*10^3 kg/m^3
		Vr = (w_cat / rho_b) * (1 - phi)  # reactor volume [m^3], calculated from w_cat.
		# Vr = 0.01e-6  # [m^3]
		### parameters end

		# read species
		species = self.get_unique_species()

		ncomp = len(species)
		ngas  = len(list(filter(lambda x: "surf" not in x, species)))

		# read pre-exponential (in [m, mol, s])
		#Afor, Afor_type = pickle.load(open("pre_exp.pickle", "rb"))
		# multipy temperature and surface density
		#for i, itype in enumerate(Afor_type):
		#	if itype=="gas":
		#		Afor[i] *= np.sqrt(T)
		#	elif itype=="ads":
		#		Afor[i] *= np.sqrt(T) / sden
		#	elif itype=="lh" or itype=="des":
		#		Afor[i] *= (T / sden)

		# temporary
		method = "emt"

		deltaE = self.get_reaction_energies(method=method)
		deltaS = self.get_entropy_differences()
		TdeltaS = T*deltaS
		TdeltaS = TdeltaS * eVtokJ

		# calculate deltaG
		deltaG = deltaE - TdeltaS

		# get Ea
		Ea = alpha * (deltaE / eVtokJ) + beta
		Ea = Ea * eVtokJ

		Kpi = np.exp(-deltaG / R / T)  # in pressure unit
		# Kci = Kpi*(101325/R/T)     # convert to concentration unit
		Kci = Kpi * (R * T / 1)  # convert to concentration unit

		tau = Vr / v0  # residence time [sec]
		# tau = 1.0  # residence time [sec]

		# output results here
		print("Ea [kJ/mol]:", Ea)
		print("deltaE [kJ/mol]:", deltaE)
		print("TdeltaS [kJ/mol]:", TdeltaS)
		print("deltaG [kJ/mol]:", deltaG)
		print("residence time [sec]: {0:5.3e}, GHSV [hr^-1]: {1:3d}".format(tau, int(60**2 / tau)))

		# now solve the ODE
		t0, tf = 0, tau
		dt = tf * 1.0e-3
		t_span = (t0, tf)
		t_eval = np.arange(t0, tf, dt)

		print("species:", species)
		# C0 = PinPa / R*T
		x0 = np.zeros(ncomp)
		x0[1] = 1.0  # CO2
		x0[2] = 4.0  # H2

		x0[-1] = 1.0  # vacancy
		C0 = Pin / (R * T * 1e3)  # density calculated from pressure. Note: R is in kJ/mol/K.

		# normalize x0 gas part
		tot = np.sum(x0[:ngas])
		for i, j in enumerate(x0):
			if i <= ngas:
				x0[i] = x0[i] / tot

		C0 *= x0

		soln = solve_ivp(fun=lambda t, C: func(t, C, rateconstants,
						Kci, T, sden, area, Vr, ngas, ncomp),
						t_span=t_span, t_eval=t_eval, y0=C0,
						rtol=1e-5, atol=1e-7, method="LSODA")  # method:BDF, Radau, or LSODA
		print(soln.nfev, "evaluations requred.")

		fig, [fig1, fig2] = plt.subplots(ncols=2, figsize=(10, 4))

		for i, isp in enumerate(species):
			if "surf" in isp:
				fig2.plot(soln.t, soln.y[i], label="theta{}".
						format(isp.replace("_", "").replace("surf", "")))
			else:
				fig1.plot(soln.t, soln.y[i], label="{}".
						format(isp))

		fig1.set_xlabel("times /s")
		fig1.set_ylabel("concentration /arb.units")
		fig2.set_xlabel("times /s")
		fig2.set_ylabel("concentration /arb.units")
		fig1.legend()
		fig2.legend()

		plt.show()

		return None

	def draw_network(self, rate=None):
		pass


class ReactionsOld:
	reactions = []

	def __init__(self, name=None):
		self.name = name
		self.file = None
		self.r_ads = None
		self.r_site = None
		self.r_coef = None
		self.p_ads = None
		self.p_site = None
		self.p_coef = None

	def count_species(self):
		pass

	def solve_ode(self):
		pass

	def add_reaction(self, rxn):
		self.reactions.append(rxn)
		print(self.reactions)

	def read_from_file(self, file):
		print("read from", file)
		self.file = file
		r_ads, r_site, r_coef, p_ads, p_site, p_coef = self._get_reac_and_prod()
		self.r_ads = r_ads
		self.r_site = r_site
		self.r_coef = r_coef
		self.p_ads = p_ads
		self.p_site = p_site
		self.p_coef = p_coef

	def read_reactionfile(self):
		lines = drop_comment_and_branck_lines(self.file)

		numlines = len(lines)

		reac = list(range(numlines))
		rxn = list(range(numlines))
		prod = list(range(numlines))

		for i, line in enumerate(lines):
			text = line.replace("\n", "").replace(">", "").split("--")
			reac_tmp = text[0]
			rxn_tmp = text[1]
			prod_tmp = text[2]

			reac[i] = re.split(" \+ ", reac_tmp)  # for cations
			prod[i] = re.split(" \+ ", prod_tmp)  # for cations

			reac[i] = remove_space(reac[i])
			prod[i] = remove_space(prod[i])

			rxn[i] = reac[i][0] + "_" + rxn_tmp

		return reac, rxn, prod

	def _get_reac_and_prod(self):
		import numpy as np
		import os
		import sys
		#
		# form reactant and product information
		#
		reac, rxn, prod = self.read_reactionfile()

		rxn_num = len(rxn)

		r_ads = list(range(rxn_num))
		r_site = [[] for i in range(rxn_num)]
		r_coef = [[] for i in range(rxn_num)]

		p_ads = list(range(rxn_num))
		p_site = list(range(rxn_num))
		p_coef = list(range(rxn_num))

		for irxn, jrnx in enumerate(rxn):
			ireac = reac[irxn]
			iprod = prod[irxn]
			ireac_num = len(ireac)
			iprod_num = len(iprod)
			#
			# reactant
			#
			r_ads[irxn] = list(range(ireac_num))
			r_site[irxn] = list(range(ireac_num))
			r_coef[irxn] = list(range(ireac_num))

			for imol, mol in enumerate(ireac):
				r_site[irxn][imol] = []
				r_ads[irxn][imol] = []
				#
				# coefficient
				#
				if "*" in mol:
					r_coef[irxn][imol] = int(mol.split("*")[0])
					rest = mol.split("*")[1]
				else:
					r_coef[irxn][imol] = 1
					rest = mol

				# site
				if ',' in rest:
					sites = rest.split(',')
					for isite, site in enumerate(sites):
						r_site[irxn][imol].append(site.split('_')[1])
						r_ads[irxn][imol].append(site.split('_')[0])
				elif '_' in rest:
					r_site[irxn][imol].append(rest.split('_')[1])
					r_ads[irxn][imol].append(rest.split('_')[0])
				else:
					r_site[irxn][imol].append('gas')
					r_ads[irxn][imol].append(rest)
			#
			# product
			#
			p_ads[irxn] = list(range(iprod_num))
			p_site[irxn] = list(range(iprod_num))
			p_coef[irxn] = list(range(iprod_num))

			for imol, mol in enumerate(iprod):
				p_site[irxn][imol] = []
				p_ads[irxn][imol] = []
				#
				# coefficient
				#
				if "*" in mol:
					p_coef[irxn][imol] = int(mol.split("*")[0])
					rest = mol.split("*")[1]
				else:
					p_coef[irxn][imol] = 1
					rest = mol

				# site
				if ',' in rest:
					sites = rest.split(',')
					for isite, site in enumerate(sites):
						p_site[irxn][imol].append(site.split('_')[1])
						p_ads[irxn][imol].append(site.split('_')[0])
				elif '_' in rest:
					p_site[irxn][imol].append(rest.split('_')[1])
					p_ads[irxn][imol].append(rest.split('_')[0])
				else:
					p_site[irxn][imol].append('gas')
					p_ads[irxn][imol].append(rest)

		# print("irxn=%d, %s-->%s, coef: %s-->%s, site:%s-->%s"
		# % (irxn, r_ads[irxn], p_ads[irxn], r_coef[irxn],
		# p_coef[irxn], r_site[irxn], p_site[irxn]))

		return r_ads, r_site, r_coef, p_ads, p_site, p_coef


def drop_comment_and_branck_lines(file):
	# drop comment and branck lines
	with open(file, "r") as f:
		lines = f.readlines()
		newlines = []
		for line in lines:
			if not (re.match(r"^#", line)) and not (re.match(r"^s*$", line)):
				newlines.append(line)

		return newlines


def remove_space(obj):
	newobj = [0] * len(obj)
	if isinstance(obj, str):
		# string
		newobj = obj.replace(" ", "")
	elif isinstance(obj, list):
		# list
		for i, obj2 in enumerate(obj):
			if isinstance(obj2, list):
				# nested list
				for ii, jj in enumerate(obj2):
					jj = jj.strip()
				newobj[i] = jj
			elif isinstance(obj2, str):
				# simple list
				obj2 = obj2.replace(" ", "")
				newobj[i] = obj2
			elif isinstance(obj2, int):
				# integer
				newobj[i] = obj2
			else:
				newobj[i] = obj2
	else:  # error
		print("remove_space: input str or list")

	return newobj
