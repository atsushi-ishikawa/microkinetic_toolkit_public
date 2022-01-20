#!/usr/bin/env python3

# formal lib
import os
import copy
import numpy as np
import pandas as pd
import pickle
import warnings
from ase.db import connect
from ase import Atoms
from typing import List
# my lib
from .reaction_for_db import ReactionForDB
from .reactions_base import ReactionsBase

R = 8.314 * 1.0e-3  # gas constant [kJ/mol/K]
eVtokJ = 96.487


class ReactionsForDB(ReactionsBase):
    """
    Set of elementary reactions.
    Overrides the method in ReactionsBase, which is empty.
    """
    @property
    def ase_db(self):
        return self._ase_db

    @ase_db.setter
    def ase_db(self, db_str: str):
        self._ase_db = connect(db_str)

    @property
    def uniq_species(self):
        if not hasattr(self, "_tot_uniq_species"):
            self._tot_uniq_species = self.get_unique_species()
        return self._tot_uniq_species

    @property
    def react_e_dict(self):
        if not hasattr(self, "_react_e_dict"):
            warning = "set react_e_dict by default state_key and select_key.\n"\
                      "state_key: chem_react_symbol, select_key: min_e_data."
            self._set_react_e_dict_from_asedb()
        return self._react_e_dict

    def set_react_e_dict_from_asedb(self,
                                    state_key="chem_react_symbol",
                                    select_key="min_e_data"):
        assert hasattr(self, "_ase_db")
        # invalid log list. 
        self._invalid_chemstr_list = []
        e_dict = {}
        self._state_key = state_key
        self._select_key = select_key
        uniq_species = self.get_unique_species()
        for uniq_sym in uniq_species:
            query_dict = {}
            query_dict[self._state_key] = str(uniq_sym)
            query_dict[self._select_key] = True
            rows_list = self._asedb.select(**query_dict),
            if len(rows_list) == 1:
                row = rows_list[0]
                E = row.energy
                if self.vibE_option:
                    zero_vE = self._get_zero_vib_e_from_row(row)
                    E = E + zero_vE
                elif G_option:
                    raise NotImplementedError("")
            elif len(rows_list) == 0:
                self._invalid_chemstr_list.append(uniq_sym)
            else:
                raise NotImplementedError("")
            e_dict[uniq_sym] = E
        self._react_e_dict = e_dict
        self.show_react_e_status()

    def show_react_e_status(self):
        if self.vibE_option:
            print("E including zero point energy.")
        elif self.G_option:
            print("E is substituted by G")
        else:
            print("E is simple Energy.")

    @property
    def vibE_option(self):
        if hasattr(self, "_vibE_option"):
            return self._vibE_option
        else:
            return False

    @property
    def G_option(self):
        if hasattr(self, "_G_option"):
            return self._G_option
        else:
            return False

    def set_vibE_option(self,
                        vibE_option: bool):
        self._vibE_option = vibE_option

    def set_G_option(self, G_option: bool,
                     T: float = None, P: float = None):
        self._G_option = G_option

    @property
    def invalid_chemstr_list(self):
        return self._invalid_chemstr_list

    def get_reaction_energies(self, zero_vib: bool = False) -> np.ndarray:
        """
        get the reaction energies (deltaEs) for all asedatabse.

        Returns:
                deltaEs: numpy array
        """
        assert hasattr(self, "_ase_db")
        deltaEs = np.zeros(len(self.reaction_list))
        for i, reaction in enumerate(self.reaction_list):
            assert isinstance(reaction, ReactionForDB)
            deltaEs[i] = reaction.get_delta_e(self.react_e_dict)
        return deltaEs

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
            reaction = ReactionForDB.from_dict(ddict)
            reaction_list.append(reaction)
        return cls(reaction_list)


class DBUtilMixin():
    def _gene_reaction_e(self) -> Iterator:
        assert hasattr(self, "_ase_db")
        for reaction in self.reaction_list:
            assert isinstance(reaction, ReactionForDB)
            deltaE = reaction.get_delta_e(self.react_e_dict)
            yield deltaE

    def _gene_react_symbols(self):
        for reaction in self.reaction_list:
            assert isinstance(reaction, ReactionForDB)
            sym = reaction._reaction_str
            yield sym
 
    def _gene_reactants_fig(self, angles: str = None):
        for reaction in self.reaction_list:
            atoms_list = self._get_atoms_reacts_from_asedb(
                                                     reaction)
    def _get_atoms_reacts_from_asedb(
                                   self,
                                   reaction: ReactionForDB) -> List[Atoms]:
        atoms_list = []
        for react_sp in reaction.reactants_species:
            query_dict = {}
            query_dict[self._state_key] = str(react_sp)
            query_dict[self._select_key] = True
            rows_list = self._asedb.select(**query_dict),
            if len(rows_list) == 1:
                row = row_list[0]
                atoms = row.toatoms()
                atoms_list.append(atoms)
            elif len(rows_list) == 0:
                warn_mes = ""
                warnings.warn()
            else:
                raise NotImplementedError("")
            return atoms_list


    def _get_atoms_products_from_asedb(
                                   self,
                                   reaction: ReactionForDB) -> List[Atoms]:
        atoms_list = []
        for react_sp in reaction.products_species:
            query_dict = {}
            query_dict[self._state_key] = str(react_sp)
            query_dict[self._select_key] = True
            rows_list = self._asedb.select(**query_dict),
            if len(rows_list) == 1:
                row = row_list[0]
                atoms = row.toatoms()
                atoms_list.append(atoms)
            elif len(rows_list) == 0:
                warn_mes = ""
                warnings.warn()
            else:
                raise NotImplementedError("")
            return atoms_list

    def _gene_products_fig(self, angles: str = None):
        raise NotImplementedError("")

    def plot_deltaE(self):
        e_iter = self._gene_reaction_e()
        sym_iter = self._gene_react_symbols()
        fig, axes = plt.subplots()
        for sym, e in zip(sym_iter, e_iter):


    def save_fig_deltaE(self, opng: str):
        raise NotImplementedError("")

    def dump_fig_deltaE(self, dump_dp: str):
        raise NotImplementedError("")


class ReactionsForDBTest(ReactionsForDB):
    def _set_react_e_dict_from_asedb(self,
                                     state_key="chem_react_symbol",
                                     select_key="min_e_data"):
        e_dict = {}
        uniq_species = self.get_unique_species()
        for uniq_sym in uniq_species:
            query_dict = {}
            query_dict[state_key] = str(uniq_sym)
            query_dict[select_key] = True
            """
            rows_list = self._asedb.select(**query_dict),
            if len(rows_list) == 1:
                row = rows_list[0]
                E = row.energy
                # zero_vE = self._get_zero_vib_e_from_row(row)
            else:
                import ipdb; ipdb.set_trace()
            """
            e_dict[uniq_sym] = 10.0
        self._react_e_dict = e_dict
