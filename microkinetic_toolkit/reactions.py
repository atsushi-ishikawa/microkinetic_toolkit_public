from microkinetic_toolkit.reaction import Reaction
import os
import numpy as np
import pandas as pd

R = 8.314 * 1.0e-3  # gas constant [kJ/mol/K]
eVtokJ = 96.487

class Reactions:
	"""
	Set of elementary reactions.
	"""
	def __init__(self, reaction_list: list):
		self.reaction_list = reaction_list
		self._ase_db = None
		self._calculator = None
		self._surface = None
		self._alpha = None
		self._beta  = None
		self._sden  = None
		self._v0    = None
		self._wcat  = None
		self._area  = None
		self._phi   = None
		self._rho_b = None
		self._Vr    = None

	def __getitem__(self, index):
		return self.reaction_list[index]

	def __len__(self):
		return len(self.reaction_list)

	@property
	def calculator(self):
		return self._calculator

	@calculator.setter
	def calculator(self, calculator_str: str):
		self._calculator = calculator_str

	@property
	def ase_db(self):
		return self._ase_db

	@ase_db.setter
	def ase_db(self, db_file: str):
		self._ase_db = db_file

	@property
	def surface(self):
		return self._surface

	@surface.setter
	def surface(self, surface):
		self._surface = surface

	def set_kinetic_parameters(self, alpha=1.0, beta=1.0, sden=1.0e-5, v0=1.0e-5, wcat=1.0e-3, phi=0.5, rho_b=1.0e3):
		"""
		Set various parameters.

		Args:
			alpha: BEP alpha
			beta: BEP beta
			sden: side density [mol/m^2]
			v0: volumetric flowrate [m^3/sec]. 1 [m^2/sec] = 1.0e6 [mL/sec] = 6.0e7 [mL/min]
			wcat: catalyst weight [kg]
			phi: porosity
			rho_b: density of catalyst [kg/m^3]. typical is 1.0 g/cm^3 = 1.0*10^3 kg/m^3
		Returns:
			None:
		"""
		self._alpha = alpha
		self._beta  = beta
		self._sden  = sden
		self._v0    = v0
		self._wcat  = wcat
		self._area  = 1000*wcat  # surface area. [m^2/kg] (e.g. BET) * [kg] --> [m^2]
		self._phi   = phi
		self._rho_b = rho_b
		self._Vr    = (wcat/rho_b)*(1-phi)  # reactor volume [m^3], calculated from w_cat.
		#self._Vr  = 0.01e-6  # [m^3]

		return None

	def to_tdb(self, db_file: str, update=False):
		tdb = TinyDB(db_file)
		for reaction in self.reaction_list:
			if update:
				reaction.update_tdb(tdb)
			else:
				reaction.to_tdb(tdb)

	def to_openfoam(self, file: str):
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

	def to_csv(self, file: str):
		"""
		Generate csv file containing elemtary reactions.

		Args:
			file: csv file name
		"""
		df = DataFrame(self._generate_reactions_dict())
		df.to_csv(file)

	@classmethod
	def from_csv(cls, csv_file: str):
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

	def get_reaction_energies(self):
		"""
		Calculate the reaction energies (deltaEs) for all the elementary reactions.

		Returns:
			deltaEs: numpy array
		"""
		deltaEs = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			deltaEs[i] = reaction.get_reaction_energy(surface=self._surface,
													  calculator=self._calculator,
													  ase_db=self._ase_db)
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
			ks[i]  = reaction.get_rate_constant(deltaE, T, alpha=self._alpha, beta=self._beta, sden=self._sden)

		return ks

	def do_microkinetics(self, deltaEs=None, ks=None, T=300.0, P=1.0, ratio=1.0):
		"""
		Do microkinetic analysis.

		Args:
			deltaEs: reaction energies.
			ks: rate constants in forward direction.
			T: temperature [K]
			P: total pressure [bar]
			ratio: pressure ratio of inlet (dict) [-]
		Returns:
			None
		"""
		if ks is None:
			raise ValueError("rate constant not found")

		odefile = "tmpode.py"
		self.make_rate_equation(odefile=odefile)
		self.solve_rate_equation(odefile=odefile, deltaEs=deltaEs, ks=ks, T=T, P=P, ratio=ratio)
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

		# template - start
		lines = [
		"\tkrev = kfor / Kc\n",
		"\ttheta = c[0:ncomp]\n",
		"\ttheta = theta * sden\n"
		]
		fout.writelines(lines)
		# template - end

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

			# making dict1, a dictionary with hash of species-number and molecule
			for side in ["reactant", "product"]:
				if side == "reactant":
					terms = reaction.reactants
					list = list_r
				elif side == "product":
					terms = reaction.products
					list = list_p
				else:
					raise ValueError("asdf")

				for term in terms:
					spe, site = term[1], term[2]
					#if site != 'gas':
					#	spe += "_surf"
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

						#if site != "gas":
						#	spe += "_surf"
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
							adsorbate = dict1[mem]
							if spe == adsorbate:
								coef = coefs[imol]

						if coef == 0:
							raise ValueError("something wrong")

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

		# tempelate - start
		lines = [
		"\tif ncomp > ngas:\n",
		"\t\trate[ngas:ncomp] = rate[ngas:ncomp]*(1/sden)  # surface\n"
		]
		# tempelate - end

		fout.writelines(lines)

		fout.write("\treturn rate\n")
		fout.write("\n")
		fout.close()

		return None

	def solve_rate_equation(self, deltaEs=None, ks=None, odefile=None, T=300.0, P=1.0, ratio=1.0):
		"""
		Solve rate equations.

		Args:
			deltaEs: reaction energies [eV]
			ks: rate constant in forward direction
			odefile: ODE file
			T: temperature [K]
			P: total pressure [bar]
			ratio: pressure ratio of inlet (dict) [-]
		Returns:
			None
		"""
		import numpy as np
		from scipy.integrate import solve_ivp
		import pickle
		from microkinetic_toolkit.tmpode import func

		if ks is None:
			raise ValueError("rate constants not found")
		if odefile is None:
			raise ValueError("ODE file not found")

		np.set_printoptions(precision=3, linewidth=100)

		Pin = P*1e5  # inlet pressure [Pascal]

		# read species
		species = self.get_unique_species()
		print("species:", species)

		ncomp = len(species)
		ngas  = len(list(filter(lambda x: "surf" not in x, species)))
		#
		# thermodynamics
		#
		deltaEs *= eVtokJ  # reaction energy
		deltaS   = self.get_entropy_differences()  # entropy
		deltaS  *= eVtokJ
		TdeltaS  = T*deltaS
		deltaG   = deltaEs - TdeltaS  # Gibbs energy

		# equilibrium constants
		Kpi = np.exp(-deltaG/R/T)  # in pressure unit
		# Kci = Kpi*(101325/R/T)   # convert to concentration unit
		Kci = Kpi*(R*T/1)          # convert to concentration unit

		# tau = self._Vr/self._v0  # residence time [sec]
		tau = 100.0  # residence time [sec]

		# empirical correction
		ks *= 1.0e10

		# output results here
		print("deltaEs [kJ/mol]:", deltaEs)
		print("TdeltaS [kJ/mol]:", TdeltaS)
		print("deltaG [kJ/mol]:", deltaG)
		print("ks [-]:", ks)
		print("res. time [sec]: {0:5.3e}, GHSV [hr^-1]: {1:3d}".format(tau, int(60**2/tau)))

		# now solve the ODE
		t0, tf = 0, tau
		dt = tf * 1.0e-3
		t_span = (t0, tf)
		t_eval = np.arange(t0, tf, dt)

		# C0 = PinPa / R*T
		x0 = np.zeros(ncomp)
		for i, j in enumerate(species):
			val = ratio.get(j)
			x0[i] = val if val is not None else 0.0

		if ncomp > ngas:
			x0[-1] = 1.0  # surface exists ... put vacancy at last

		# normalize x0 gas part
		tot = np.sum(x0[:ngas])
		for i, j in enumerate(x0):
			if i <= ngas:
				x0[i] = x0[i] / tot

		C0 = Pin/(R*T*1e3)  # density calculated from pressure. Note: R is in kJ/mol/K.
		C0 *= x0

		soln = solve_ivp(fun=lambda t, C: func(t, C, ks,
						Kci, T, self._sden, self._area, self._Vr, ngas, ncomp),
						t_span=t_span, t_eval=t_eval, y0=C0,
						rtol=1e-5, atol=1e-7, method="LSODA")  # method:BDF, Radau, or LSODA
		print(soln.nfev, "evaluations requred.")

		self.draw_molar_fraction_change(soln=soln, showfigure=True, savefigure=False)
		return None

	def draw_molar_fraction_change(self, soln=None, showfigure=False, savefigure=False, filename="result.png"):
		"""
		Draw molar fraction change with time.

		Args:
			soln: solution from solve_ivp
			showfigure: whether to show figure
			savefigure: whether to save figure
			filename: file name when saving figure
		"""
		import matplotlib.pyplot as plt

		if soln is None:
			raise Exception("Nothing to plot")

		species = self.get_unique_species()

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

		if showfigure:
			plt.show()
		if savefigure:
			plt.savefig(filename)

		return None

	def draw_network(self, rate=None):
		pass

	def do_preparation(self):
		from .preparation import preparation
		# should pass the adsorbate set
		preparation.prepare(self.get_unique_species())
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

	def read_from_file(self, file):
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
