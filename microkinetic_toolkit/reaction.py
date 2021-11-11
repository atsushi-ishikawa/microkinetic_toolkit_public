import os
import sys
import re
import textwrap
import numpy as np
import pandas as pd
from pandas import DataFrame
from tinydb import TinyDB, Query

class Reaction:
	"""
	Elemetary reaction class.
	"""
	class IdIterator(object):
		def __init__(self):
			self.num = 0

		def __next__(self):
			num = self.num
			self.num += 1
			return num

	id_iterator = IdIterator()

	def __init__(self, reaction_str, reaction_id=None):
		self._reaction_str = reaction_str
		if reaction_id is None:
			reaction_id = next(self.id_iterator)
		self._reaction_id = reaction_id
		self._parse_reaction_str()

	def _parse_reaction_str(self):
		reactants_str, products_str = self._reaction_str.split("->")
		self.reactants = self._from_halfside_to_nsp_list(reactants_str)
		self.products  = self._from_halfside_to_nsp_list(products_str)

	def _from_halfside_to_nsp_list(self, halfside_str):
		terms = halfside_str.split("+")
		list_nspecies = [self._change_term_to_coef_and_species(term) for term in terms]
		return list_nspecies

	def _change_term_to_coef_and_species(self, term):
		term = term.strip()

		# coefficient
		re_coef = re.compile("^[0-9]+")
		re_coef_search = re_coef.search(term)

		# site
		re_site = re.compile("_(atop|fcc|hcp)")
		re_site_search = re_site.search(term)

		# coef
		if re_coef_search is None:  # has no coefficient
			coef = 1
		else:  # has coefficient
			coef = int(re_coef_search[0])

		# site
		if re_site_search is None:  # has no site -- gas
			site = "gas"
		else:  # has site
			site = re_coef_search[1]

		spe = re_site.sub("", re_coef.sub("", term))
		spe = spe.strip()

		return (coef, spe, site)

	def to_dict(self):
		"""
		convert to dict.
		Returns: base_attributed dict
		"""
		def extract_attribute_dict(instance, keys, default=None):
			ddict = {}
			for key in keys:
				val = getattr(instance, key, default)
				ddict[key] = val
			return ddict
		return extract_attribute_dict(self, self.base_attributed_keys)

	@classmethod
	def from_dict(cls, ddict):
		reaction_str = ddict["_reaction_str"]
		reaction_id  = ddict["_reaction_id"]
		return cls(reaction_str, reaction_id)

	# tinyDB
	def to_tdb(self, tdb):
		tdb.insert(self.to_dict())

	def update_tdb(self, tdb):
		ddict = self.to_dict()
		query = Query()
		q_ins = getattr(query, "_reaction_id")
		tdb.update(ddict, q_ins == ddict["_reaction_id"])

	@property
	def base_attributed_keys(self):
		base_keys = ["_reaction_str", "_reaction_id"]
		return base_keys

	@property
	def unique_species(self):
		"""
		Get unique species in an elementary reaction.
		Returns:
			Unique chemical species -- list of strings.
		"""
		reac = self.reactants_species
		prod = self.products_species
		return list(set(reac).union(set(prod)))

	@property
	def reactants_species(self):
		"""
		Get reactant chemical species.
		Returns:
			reactant species in the list of strings
		"""
		lst = []
		for i in self.reactants:
			lst.append(i[1])
		return lst

	@property
	def products_species(self):
		"""
		Get product chemical species.
		Returns:
			product species in the list of strings
		"""
		lst = []
		for i in self.products:
			lst.append(i[1])
		return lst

	# openfoam
	def _get_openfoam_react_str(self):
		str = {"reac": None, "prod": None}
		for side in ["reac", "prod"]:
			from_ = self.reactants if side == "reac" else self.products
			terms = []
			for coef, species in from_:
				term = species if coef == 1 else "{}{}".format(coef, species)
				terms.append(term)
			str[side] = " + ".join(terms)
		openfoam_reaction_str = str["reac"] + " = " + str["prod"]
		return openfoam_reaction_str

	def to_openfoam_paramdict(self):
		"""
		Get openfoam paramter as dict.
		Returns:
			dict
		"""
		param_dict = {}
		param_dict["type"] = "TestReactionType"
		param_dict["reaction"] = self._get_openfoam_react_str()
		param_dict["A"] = 1.0e16
		param_dict["beta"] = 0
		param_dict["Ta"] = 10000
		return param_dict

	def adsorbate_on_surface(self, surface=None, height=0.0):
		"""
		Adsorbate molecule on surface.
		Args:
			surface:
			height:
		Returns:
			atoms_dict: {"reactants": ASE Atoms for reactants, "products": for products}
		"""
		from ase.build import fcc111, add_adsorbate, molecule
		from ase.visualize import view

		atoms_dict = {"reactants": [], "products": []}

		if surface is None:
			surface = fcc111("Au", size=[2, 2, 3], vacuum=10.0)
		if height == 0.0:
			height = 2.0

		for side in ["reactants", "products"]:
			sequence  = self.reactants if side == "reac" else self.products
			for i in sequence:
				_, mol, _ = i
				ads = molecule(mol)
				surf_copy = surface.copy()

				# adjust height
				shift = min(ads.positions[:, 2])
				height -= shift
				add_adsorbate(surf_copy, ads, offset=(0, 0), position=(0, 0), height=height)
				surf_copy.pbc = True
				#view(surf_copy)

				atoms_dict[side].append(surf_copy)

		return atoms_dict

	def get_reaction_energy(self, surface=None, calculator=None, ase_db=None):
		"""
		Calculate reaction energy for an elementary reaction.
		Args:
			surface: Atoms
			calculator: currently Vasp or EMT
			ase_db: ASE database file for existing calculation
			        Expects energy is stored by
			        db.write(atoms, data={"energy": -123.4, ...})
		Returns:
			deltaE (float)
		"""
		from ase.calculators.emt import EMT

		def search_energy_from_ase_db(atom, ase_db):
			from ase.db import connect

			db = connect(ase_db)
			formula = atom.get_chemical_formula()
			try:
				energy = db.get(name=formula).data.energy
			except:
				print("{0:} is not found in {1:}.".format(formula, db))
				energy = 0.0

			return energy

		# set calculator
		if calculator == "vasp":
			# TODO: do vasp setting
			pass
		elif calculator == "emt":
			calc = EMT()
		else:
			raise Exception("use vasp or emt for calculator")

		atoms_dict = self.adsorbate_on_surface(surface=surface)
		energy_dict = {"reactants": 0.0, "products": 0.0}

		for side in ["reactants", "products"]:
			for iatom in atoms_dict[side]:
				if ase_db is not None:
					energy = search_energy_from_ase_db(iatom, ase_db)
				else:
					iatom.calc = calc
					energy = iatom.get_potential_energy()

				energy_dict[side] += energy

		deltaE = energy_dict["products"] - energy_dict["reactants"]

		return deltaE

	def get_entropy_difference(self):
		"""
		Calculate entropy difference (deltaS, in eV/K) along the reaction.
		Returns:
			deltaS (float)
		"""
		from ase.collections import methane

		reac_S = 0.0
		prod_S = 0.0

		# reactants
		for mol in self.reactants:
			spe, site = mol[0], mol[1]

			if site != 'gas' or spe == 'surf' or spe == 'def':
				# surface species
				entropy = 0.0
			else:
				try:
					entropy = methane.data[spe]['molecular_entropy']
				except:
					entropy = 1.0e-3  # rough entropy estimation ... 1 meV/K

			reac_S += entropy

		# products
		for mol in self.products:
			spe, site = mol[0], mol[1]

			if site != 'gas' or spe == 'surf' or spe == 'def':
				entropy = 0.0
			else:
				try:
					entropy = methane.data[spe]['molecular_entropy']
				except:
					entropy = 0.0

			prod_S += entropy

		deltaS = np.sum(prod_S) - np.sum(reac_S)
		return deltaS

	def get_rate_constant(self, deltaE=0.0, T=300.0):
		"""
		Calculate rate constant from reaction energy (deltaE).
		Args:
			deltaE: reaction energy [eV]
			T: temperature [K]
		Returns:
			rate constant
		"""
		import numpy as np
		from ase.units import create_units

		units = create_units('2014')
		amu = units['_amu']  # 1.66e-27 [kg]
		kB = units['_k']  # 1.38e-23 [J/K]
		hplanck = units['_hplanck']  # 6.63e-34 [J*s]
		Nav = units['_Nav']

		type = "LH"
		sden = 1.0e10  # site density

		# calculate Ea from deltaE
		alpha = 1.0
		beta  = 1.0
		Ea = alpha*deltaE + beta

		if type == "LH":
			A = kB*T/hplanck/sden  # [s^-1]

		exp = np.exp(-Ea/kB*T)

		rateconst = A*exp

		return rateconst
