import re
import textwrap
import pandas as pd
from pandas import DataFrame
from tinydb import TinyDB, Query
from sympy import Symbol

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
		self.products = self._from_hside_to_nsp_list(products_str)

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

def extract_attribute_dict(instance, keys, default=None):
	ddict = {}
	for key in keys:
		val = getattr(instance, key, default)
		ddict[key] = val
	return ddict


class Reactions:
	def __init__(self, react_data_list):
		self.react_data_list = react_data_list

	def __getitem__(self, index):
		return self.react_data_list[index]

	def to_tdb(self, db_file, update=False):
		tdb = TinyDB(db_file)
		for react_data in self.react_data_list:
			if update:
				react_data.update_tdb(tdb)
			else:
				react_data.to_tdb(tdb)

	def to_openfoam(self, file):
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
		species_set = set([])
		for react_data in self.react_data_list:
			species_set.update(react_data.unique_species)
		return species_set

	def _get_openfoam_reactions_info(self):
		ini = "reaction\n{\n"
		yield ini
		for react_info in self._get_each_react_info():
			react_info = textwrap.indent(react_info, prefix="   ")
			yield react_info
		fin = "}\n"
		yield fin

	def _get_each_react_info(self):
		for react_data in self.react_data_list:
			paramdict = react_data.to_openfoam_paramdict()
			ini = "TestReaction{}\n".format(react_data._reaction_id)
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
		for react_data in self.react_data_list:
			ddict = react_data.to_dict()
			yield ddict

	def to_csv(self, file):
		df = DataFrame(self._generate_reactions_dict())
		df.to_csv(file)

	@classmethod
	def from_csv(cls, csv_file):
		df = pd.read_csv(csv_file, index_col=0)
		react_data_list = []
		for i, row in df.iterrows():
			ddict = row.to_dict()
			react_data = Reaction.from_dict(ddict)
			react_data_list.append(react_data)
		return cls(react_data_list)

class ReactionsOld:
	"""
	Set of elementary reactions.
	"""
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
