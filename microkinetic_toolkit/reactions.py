from microkinetic_toolkit.reaction import Reaction
import os
import copy
import numpy as np
import pandas as pd
import pickle

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
		self._bep_param = None
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
		# shift
		zmin = np.min(surface.positions[:, 2])
		surface.translate([0, 0, -zmin + 0.1])

		# sort
		surface = self.sort_atoms_by_z(surface)
		self._surface = surface

	def set_kinetic_parameters(self, bep_param=None, sden=1.0e-5, v0=1.0e-5, wcat=1.0e-3, phi=0.5, rho_b=1.0e3):
		"""
		Set various parameters.

		Args:
			bep_param: BEP alpha and beta (dict)
			sden: side density [mol/m^2]
			v0: volumetric flowrate [m^3/sec]. 1 [m^2/sec] = 1.0e6 [mL/sec] = 6.0e7 [mL/min]
			wcat: catalyst weight [kg]
			phi: porosity
			rho_b: density of catalyst [kg/m^3]. typical is 1.0 g/cm^3 = 1.0*10^3 kg/m^3
		Returns:
			None:
		"""
		if bep_param is None:
			self._bep_param = {"alpha", 1.0, "beta", 1.0}
		else:
			self._bep_param = bep_param
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
		species_set = list(species_set)
		species_set = sorted(species_set)

		# move "surf" to last
		if "surf" in species_set:
			species_set.remove("surf")
			species_set.append("surf")

		return list(species_set)

	def get_index_of_species(self, specie):
		"""
		Get index of species.

		Args:
			specie: string
		Returns:
			specie_index: int
		"""
		species = self.get_unique_species()
		specie_index = species.index(specie)
		return specie_index

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

	def freeze_surface(self):
		def set_tags(surface):
			tags = surface.get_tags()
			surface_copy = copy.deepcopy(surface)
			maxval = max(tags)
			newtags = list(map(lambda x: 0 if x > maxval//2 else 1, tags))
			newtags = np.array(newtags)
			surface_copy.set_tags(newtags)
			return surface_copy

		import ase.constraints
		surface = copy.deepcopy(self._surface)
		surface = set_tags(surface)
		c = ase.constraints.FixAtoms(indices=[atom.index for atom in surface if atom.tag == 0])
		surface.set_constraint(c)
		return surface

	def get_reaction_energies(self):
		"""
		Calculate the reaction energies (deltaEs) for all the elementary reactions.

		Returns:
			deltaEs: numpy array
		"""
		# freeze surface here
		surface = self.freeze_surface()

		deltaEs = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			deltaEs[i] = reaction.get_reaction_energy(surface=surface, calculator=self._calculator,
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
			kfor: forward rate constants (numpy array)
		"""
		kfor = np.zeros(len(self.reaction_list))
		for i, reaction in enumerate(self.reaction_list):
			index  = reaction._reaction_id
			deltaE = deltaEs[index]
			kfor[i]  = reaction.get_rate_constant(deltaE, T, bep_param=self._bep_param, sden=self._sden)

		return kfor

	def do_microkinetics(self, deltaEs=None, kfor=None, T=300.0, P=1.0, ratio=1.0):
		"""
		Do microkinetic analysis.

		Args:
			deltaEs: reaction energies.
			kfor: rate constants in forward direction.
			T: temperature [K]
			P: total pressure [bar]
			ratio: pressure ratio of inlet (dict) [-]
		Returns:
			None
		"""
		if kfor is None:
			raise ValueError("rate constant not found")

		odefile = "tmpode.py"
		self.make_rate_equation(odefile=odefile)
		self.solve_rate_equation(odefile=odefile, deltaEs=deltaEs, kfor=kfor, T=T, P=P, ratio=ratio)
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

		for irxn, reaction in enumerate(self.reaction_list):
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
					#spe_num = self.get_unique_species().index(spe)
					spe_num = self.get_index_of_species(spe)
					list.append(spe_num)
					dict1[spe_num] = spe
			# done for dict1

			# making dict2
			for side in ["reactant", "product"]:
				for direction in ["forward", "reverse"]:
					if side == "reactant" and direction == "forward":
						mol_list1 = reaction.reactants
						current_list = list_r
						mol_list2 = reaction.reactants  # list corresponding to current_list
						coefs = [i[0] for i in mol_list2]
						term = "kfor[" + str(irxn) + "]"
						sign = " - "
					elif side == "reactant" and direction == "reverse":
						mol_list1 = reaction.products
						current_list = list_r
						mol_list2 = reaction.reactants
						coefs = [i[0] for i in mol_list2]
						term = "krev[" + str(irxn) + "]"
						sign = " + "
					elif side == "product" and direction == "forward":
						mol_list1 = reaction.reactants
						current_list = list_p
						mol_list2 = reaction.products
						coefs = [i[0] for i in mol_list2]
						term = "kfor[" + str(irxn) + "]"
						sign = " + "
					elif side == "product" and direction == "reverse":
						mol_list1 = reaction.products
						current_list = list_p
						mol_list2 = reaction.products
						coefs = [i[0] for i in mol_list2]
						term = "krev[" + str(irxn) + "]"
						sign = " - "

					# making single term
					for mol in mol_list1:
						coef, spe, site = mol[0], mol[1], mol[2]
						spe_num = self.get_unique_species().index(spe)

						if site == "gas" and spe != "surf":
							theta = "c[" + str(spe_num) + "]"      # gas-phase molecule
						else:
							theta = "theta[" + str(spe_num) + "]"  # adsorbed species or bare surface

						power = coef
						if power != 1:
							theta += "**" + str(power)

						term += "*" + theta

					for i in current_list:
						if dict1[i] == "surf":
							continue  # bare surface ... skip

						coef = 0
						for imol, mol in enumerate(mol_list2):
							spe = mol[1]
							adsorbate = dict1[i]
							if spe == adsorbate:
								coef = coefs[imol]

						if coef == 0:
							raise ValueError("something wrong")

						sto_coef = str(float(coef))

						if i in dict2:
							dict2[i] += sign + sto_coef + "*" + term  # NEGATIVE
						else:
							dict2[i] = sign + sto_coef + "*" + term

						if "theta" in dict2[i]:
							dict2[i] += "*" + "(area/Vr)"
					# loop for current_list
				# direction
			# side
		# loop over elementary reactions

		# term for the vacant site
		#   formation/consumption rates for all the surface species should be
		#   subtracted from the vacant site.
		if "surf" in dict1.values():
			# surface reaction
			tmp = ""
			ind = [key for key, value in dict1.items() if value == "surf"][0]  # index for "surf"

			for imol, mol in enumerate(dict1):
				comp = dict1[imol]
				if "surf" in comp and comp != "surf":  # surface reaction but not surface itself
					tmp += " -rate[" + str(imol) + "]"

			dict2[ind] = tmp

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

	def get_equilibrium_constant_from_deltaEs(self, deltaEs=None, T=300.0):
		# get deltaG
		deltaEs  = deltaEs.copy()*eVtokJ
		deltaSs  = self.get_entropy_differences()
		deltaSs  = deltaSs.copy()*eVtokJ
		TdeltaSs = T*deltaSs
		deltaGs = deltaEs - TdeltaSs

		self.print_for_all_the_reactions(deltaGs)

		# equilibrium constants
		Kp = np.exp(-deltaGs/R/T)  # in pressure unit
		#Kc = Kp*(101325/R/T)    # convert to concentration unit (when using atm)
		Kc = Kp*((R*1.0e3)*T/1)  # convert to concentration unit (note: R is in kJ/mol/K)
		return Kc

	def solve_rate_equation(self, deltaEs=None, kfor=None, odefile=None, T=300.0, P=1.0, ratio=1.0):
		"""
		Solve rate equations.

		Args:
			deltaEs: reaction energies [eV]
			kfor: rate constant in forward direction
			odefile: ODE file
			T: temperature [K]
			P: total pressure [bar]
			ratio: pressure ratio of inlet (dict) [-]
		Returns:
			None
		"""
		import numpy as np
		from scipy.integrate import solve_ivp
		from microkinetic_toolkit.tmpode import func

		if kfor is None:
			raise ValueError("rate constants not found")
		if odefile is None:
			raise ValueError("ODE file not found")

		np.set_printoptions(precision=3, linewidth=100)

		Pin = P*1.0e5  # inlet pressure converted from bar to Pascal

		# read species
		species = self.get_unique_species()
		print("species:", species)

		ncomp = len(species)
		ngas  = len(list(filter(lambda x: "surf" not in x, species)))

		Kc = self.get_equilibrium_constant_from_deltaEs(deltaEs, T=T)

		tau = self._Vr/self._v0  # residence time [sec]
		# tau = 1.0  # residence time [sec]

		# empirical correction
		#kfor *= 1.0e0

		# output results here
		print("Kc [-]:", Kc)
		print("kfor [-]:", kfor)
		print("residence time [sec]: {0:5.3e}, GHSV [hr^-1]: {1:3d}".format(tau, int(60**2/tau)))

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

		# method:BDF, Radau, or LSODA
		soln = solve_ivp(fun=lambda t, C: func(t, C, kfor, Kc, T, self._sden, self._area, self._Vr, ngas, ncomp),
						 t_span=t_span, t_eval=t_eval, y0=C0, rtol=1e-5, atol=1e-7, method="LSODA")
		print(soln.nfev, "evaluations requred.")

		self.draw_molar_fraction_change(soln=soln, showfigure=True, savefigure=False)
		self.save_coverage(soln=soln)
		return None

	def save_coverage(self, soln=None):
		"""
		Save surface coverage: for time-independent, output the coverage at the last time step

		Args:
			soln: solution of ivp solver
		Returns:
		"""
		import h5py

		# time dependent coverage for graph
		#dT = 10  # one record in dT points
		fac = 1.0e-6

		species = self.get_unique_species()
		ncomp = len(species)
		ngas  = len(list(filter(lambda x: "surf" not in x, species)))
		maxtime = len(soln.t)

		tcov = []  # time-dependent coverage: tcov[species][time]
		for i in range(ngas):
			tcov.append(soln.y[i])
		for i in range(ngas, ncomp):
			tcov.append(soln.y[i]*fac)  # multiply scaling factor for surface species
		tcov = np.array(tcov)

		h5file = h5py.File("coverage.h5", "w")
		h5file.create_dataset("time", shape=(maxtime,), dtype=np.float, data=soln.t)
		h5file.create_dataset("concentration", (ncomp, maxtime), dtype=np.float, data=tcov)
		h5file.close()

		#f = open(coveragefile, "wt")
		#f.write("      name      num     conc        time\n")  # header

		## take averaged value
		#with h5py.File("coverage.h5", "a") as file:
		#	for isp in range(ncomp):
		#		for itime in range(maxtime):
		#			file["concentration"][isp,itime] = tcov[isp][itime]
		#	#f.write("{0:16.14s}{1:03d}{2:12.4e}{3:12.4e}\n".format(species[isp],isp,tcov[isp][it],soln.t[it]))

		## add RXN nodes at last
		#for i in range(Nrxn):
		#	f.write("R{0:>03d}            {1:03d}{2:12.4e}{3:12.4e}\n".format(i, i, 1.0e-20, 0.0))
		#f.close()
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

	def get_rate_for_graph(self, T=300.0):
		import h5py

		conc_h5  = "coverage.h5"
		deltaEs_pickle = "deltaEs.pickle"

		rxn_num = len(self.reaction_list)

		h5file = h5py.File(conc_h5, "r")
		time = h5file["time"][:]
		cov = h5file["concentration"]

		#for iline, line in enumerate(lines):
		#	name = str(line.split()[0])
		#	if "R" not in name:
		#		spe = int(line.split()[1]) - 1  # Matlab --> python
		#		conc = float(line.split()[2])
		#		time = float(line.split()[3])
		#
		#		cov[spe].append([name, conc, time])
 		#	else:
 		#		Nrxn += 1

		max_time = len(time)

		# read deltaEs
		with open(deltaEs_pickle, "rb") as file:
			deltaEs = pickle.load(file)

		kfor = self.get_rate_constants(deltaEs=deltaEs, T=T)
		Kc = self.get_equilibrium_constant_from_deltaEs(deltaEs=deltaEs)
		krev = kfor / Kc

		AoverV = self._area / self._Vr

		# output information for graph
		rates = {"forward": None, "reverse": None}
		for direction in ["forward", "reverse"]:
			rate = np.zeros((max_time, rxn_num))
			for irxn, reaction in enumerate(self.reaction_list):
				if direction == "forward":
					sequence = reaction.reactants
					k = kfor
				elif direction == "reverse":
					sequence = reaction.products
					k = krev

				for istep in range(max_time):
					conc_all = 1.0
					for mol in sequence:
						coef, spe, site = mol

						spe_num = self.get_index_of_species(spe)
						conc = cov[spe_num][istep]
						if site != "gas":  # surface species
							conc *= self._sden * AoverV
						if coef != 1:
							conc = conc ** coef

						conc_all = conc_all * conc

					rate[istep][irxn] = k[irxn] * conc_all

			rates[direction] = rate

		# time-dependent case
		#edgefile = "edges.txt"
		#f = open(edgefile, "w")
		#f.write("    from       to         rate         time\n")

		#for istep in range(max_time):
		#	for irxn, reaction in enumerate(self.reaction_list):
		#		for direction in ["forward", "reverse"]:
		#			if direction == "forward":
		#				sequence = reaction.reactants
		#			elif direction == "reverse":
		#				sequence = reaction.products
#
#					now = 'R%03d' % irxn  # current reaction
#
#					for mol in sequence:
#						spe = mol[1]
#						if direction == "forward":
#							f.write("%10s%10s%12.4e%12.4e\n" % (spe, now, rates[direction][istep][irxn], time[istep]))
#						elif direction == "product":
#							f.write("%10s%10s%12.4e%12.4e\n" % (now, spe, rates[direction][istep][irxn], time[istep]))
#
#		f.close()
#
 		# time-independent case - write last rate
		total_rate = np.zeros(rxn_num)
		for irxn in range(rxn_num):
			total_rate[irxn] = rates["forward"][-1][irxn] - rates["reverse"][-1][irxn]

		self.print_for_all_the_reactions(total_rate)

		with open("rate.pickle", "wb") as file:
			pickle.dump(total_rate, file)

	def print_for_all_the_reactions(self, property=None):
		print("{0:^40.35s}{1:^12s}".format("reaction", "value"))
		for irxn, reaction in enumerate(self.reaction_list):
			print("{0:<40.35s}{1:>10.2e}".format(reaction._reaction_str, property[irxn]))

	def draw_network(self):
		import matplotlib.pyplot as plt
		import networkx as nx
		from math import log10
		import h5py

		label_rxn = True
		directed = True

		eps        = 1.0e-10
		edge_scale = 0.2   # scaling for edge thickness
		edge_min   = 1.0   # minimal size of edge
		rate_thre  = -15.0  # Threshold for rxn rate in log scale. Rxn with smallter than this value is discarded.

		node_scale = 100.0
		node_min   = 10.0
		surf_scale = 1.0

		rxn_num = len(self.reaction_list)

		with open("rate.pickle", "rb") as file:
			rate = pickle.load(file)

		value = np.zeros(rxn_num)
		for irxn, reaction in enumerate(self.reaction_list):
			if abs(rate[irxn]) >= rate_thre:
				value[irxn] = edge_scale * log10(abs(rate[irxn])) + edge_min
			else:
				value[irxn] = 1.0

		c_siz = 200
		c_col = "blue"
		r_siz = 10
		r_col = "black"

		coverage = True
		conc_h5 = "coverage.h5"
		if coverage:
			h5file = h5py.File(conc_h5, "r")
			cov = h5file["concentration"]
		cov = cov[:, -1]  # last

		if directed:
			G = nx.DiGraph()
		else:
			G = nx.Graph()

		for irxn, reaction in enumerate(self.reaction_list):
			rxn = reaction._reaction_str
			G.add_node(rxn, size=r_siz, color=r_col, typ='rxn')
			for direction in ["forward", "reverse"]:
				sequence = reaction.reactants if direction == "forward" else reaction.products
				for mol in sequence:
					if coverage:
						spe = mol[1]
						spe_num = self.get_index_of_species(spe)
						size = cov[spe_num] if cov[spe_num] > eps else eps
						size = node_scale * (node_min + log10(size))
						if "surf" in spe:
							size = size * surf_scale
						size = int(size)
					else:
						size = c_siz

					G.add_node(spe, size=size, color=c_col, typ='comp')

					if directed:
						if direction == "forward":
							if value[irxn] > 0:
								G.add_edge(spe, rxn, weight=abs(value[irxn]))
							else:
								G.add_edge(rxn, spe, weight=abs(value[irxn]))
						else:  # reverse
							if value[irxn] > 0:
								G.add_edge(rxn, spe, weight=abs(value[irxn]))
							else:
								G.add_edge(spe, rxn, weight=abs(value[irxn]))
					else:  # non-directed
						#G.add_edge(reac[i][ireac], rxn, weight=abs(value[i]))
						G.add_edge(spe, rxn, weight=abs(value[irxn]))

		# drawing
		siz = nx.get_node_attributes(G, 'size')
		col = nx.get_node_attributes(G, 'color')

		pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='fdp')  # prog='neato' is also a good choice
		nx.draw_networkx_nodes(G, pos, nodelist=list(siz.keys()),
							   node_size=list(siz.values()), node_color=list(col.values()), alpha=0.5)
		edges = G.edges()
		weights = [G[u][v]['weight'] for u, v in edges]

		nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.8, width=weights)

		# compound labels
		if directed:
			Gcomp = nx.DiGraph()
		else:
			Gcomp = nx.Graph()

		for n, typ in G.nodes.data('typ'):
			if typ == 'comp':
				Gcomp.add_node(n)

		nx.draw_networkx_labels(Gcomp, pos, font_size=12, font_family='Gill Sans MT')

		# reaction labels
		if label_rxn:
			rxn_label_size = 10
		else:
			rxn_label_size = 0

		if directed:
			Grxn = nx.DiGraph()
		else:
			Grxn = nx.Graph()
		#
		for n, typ in G.nodes.data('typ'):
			if typ == 'rxn':
				Grxn.add_node(n)

		nx.draw_networkx_labels(Grxn, pos, font_size=rxn_label_size, font_family='Gill Sans MT')

		plt.xticks([])
		plt.yticks([])
		plt.show()

		plt.figure(figsize=(16, 10))
		# plt.savefig("graph.eps", format="eps")

		nx.write_gexf(G, 'test.gexf')

	def do_preparation(self):
		from .preparation import preparation
		# should pass the adsorbate set
		preparation.prepare(self.get_unique_species())
		pass

	def sort_atoms_by_z(self, atoms):
		import ase
		import copy

		dtype = [("idx", int), ("z", float)]

		# get set of chemical symbols
		atoms_copy = copy.deepcopy(atoms)
		symbols = atoms_copy.get_chemical_symbols()
		elements = sorted(set(symbols), key=symbols.index)
		num_elem = []
		for i in elements:
			num_elem.append(symbols.count(i))

		iatm = 0
		newatoms = ase.Atoms()
		for inum in num_elem:
			zlist = np.array([], dtype=dtype)
			for idx in range(inum):
				tmp = np.array([(iatm, atoms_copy[iatm].z)], dtype=dtype)
				zlist = np.append(zlist, tmp)
				iatm += 1

			zlist = np.sort(zlist, order="z")

			for i in zlist:
				idx = i[0]
				newatoms.append(atoms[idx])

		newatoms.tags = atoms_copy.get_tags()
		newatoms.pbc = atoms_copy.get_pbc()
		newatoms.cell = atoms_copy.get_cell()

		return newatoms


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
		r_site = [[] for _ in range(rxn_num)]
		r_coef = [[] for _ in range(rxn_num)]

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
