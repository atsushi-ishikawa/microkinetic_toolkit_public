import os
import sys
import re
import textwrap
import pandas as pd
from pandas import DataFrame
from tinydb import TinyDB, Query
from sympy import Symbol

def extract_attribute_dict(instance, keys, default=None):
	ddict = {}
	for key in keys:
		val = getattr(instance, key, default)
		ddict[key] = val
	return ddict

class Reaction:
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
		self.sympy_reactants = self.reactants
		self.sympy_products  = self.products

	def _parse_reaction_str(self):
		reactants_str, products_str = self._reaction_str.split("->")
		self.reactants = self._from_hside_to_nsp_list(reactants_str)
		self.products  = self._from_hside_to_nsp_list(products_str)

	def _from_hside_to_nsp_list(self, hside_str):
		terms = hside_str.split("+")
		list_nspecies = [self._change_term_to_coef_and_species(term) for term in terms]
		return list_nspecies

	def _change_term_to_coef_and_species(self, term):
		term = term.strip()
		reg = re.compile("^[0-9]+")
		re_match = reg.match(term)
		if re_match is None:  # has no coefficient
			sp = term.strip()
			return (1, sp)
		else:  # has coefficient
			coef = int(re_match[0])
			sp = reg.sub("", term)
			sp = sp.strip()
			return coef, sp

	def to_dict(self):
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

	def _get_sympy_hside(self, list_nspecies):
		hside_expr = 0.0
		for num, sym in list_nspecies:
			if " *" in sym:
				sym = sym.repace("*", "{a}")
			species = Symbol(sym)
			term = num * species
			hside_expr += term
		return hside_expr

	@property
	def sympy_reactants(self):
		return self._sympy_reactants

	@sympy_reactants.setter
	def sympy_reactants(self, reactants):
		self._sympy_reactants = self._get_sympy_hside(reactants)

	@property
	def sympy_products(self):
		return self._sympy_products

	@sympy_products.setter
	def sympy_products(self, products):
		self._sympy_products = self._get_sympy_hside(products)

	@property
	def unique_species(self):
		return self.reactants_species.union(self.products_species)

	@property
	def reactants_species(self):
		return self.sympy_reactants.free_symbols

	@property
	def products_species(self):
		return self.sympy_products.free_symbols

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
		param_dict = {}
		param_dict["type"] = "TestReactionType"
		param_dict["reaction"] = self._get_openfoam_react_str()
		param_dict["A"] = 1.0e16
		param_dict["beta"] = 0
		param_dict["Ta"] = 10000
		return param_dict

	def adsorbate_on_surface(self, surface=None, height=0.0):
		from ase.build import fcc111, add_adsorbate, molecule
		from ase.visualize import view

		atoms = {"reactants": [], "products": []}

		if surface is None:
			surface = fcc111("Au", size=[2, 2, 3], vacuum=10.0)
		if height == 0.0:
			height = 2.0

		for side in ["reactants", "products"]:
			sequence  = self.reactants if side == "reac" else self.products
			for i in sequence:
				_, mol = i
				ads = molecule(mol)
				surf_copy = surface.copy()

				# adjust height
				shift = min(ads.positions[:, 2])
				height -= shift
				add_adsorbate(surf_copy, ads, offset=(0, 0), position=(0, 0), height=height)
				surf_copy.pbc = True
				#view(surf_copy)

				atoms[side].append(surf_copy)

		return atoms

	def get_reaction_energy(self, surface=None, method=None):
		from ase.calculators.emt import EMT
		calc = EMT()

		atoms = self.adsorbate_on_surface(surface=surface)
		# reactant
		r_energy = 0.0
		for i in atoms["reactants"]:
			i.calc = calc
			e = i.get_potential_energy()
			r_energy += e

		# product
		p_energy = 0.0
		for i in atoms["products"]:
			i.calc = calc
			e = i.get_potential_energy()
			p_energy += e

		deltaE = p_energy - r_energy

		return deltaE

	def get_rate_constant(self, deltaE=0.0, T=300.0):
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
