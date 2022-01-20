#!/usr/bin/env python3

# formal lib
import os
import copy
import numpy as np
import pandas as pd
import pickle
from ase.db import connect
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
            self._set_react_e_dict_from_asedb()
        return self._react_e_dict

    def _set_react_e_dict_from_asedb_v2(self,
                                        state_key="chem_react_symbol",
                                        select_key="min_e_data"):
        assert hasattr(self, "_ase_db")
        e_dict = {}
        uniq_species = self.get_unique_species()
        for uniq_sym in uniq_species:
            query_dict = {}
            query_dict[state_key] = str(uniq_sym)
            query_dict[select_key] = True
            rows_list = self._asedb.select(**query_dict),
            if len(rows_list) == 1:
                row = rows_list[0]
                E = row.energy
                # zero_vE = self._get_zero_vib_e_from_row(row)
            elif len(rows_list) == 0:
                raise NotImplementedError("")
            else:
                import ipdb
                ipdb.set_trace()
            e_dict[uniq_sym] = E
        self._react_e_dict = e_dict

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

    def _get_zero_vib_e_from_row(self, row):
        vib_e_key = "vib_zero_e"
        ex_dict = row.data
        if zero_vib_e_key in ex_dict:
            zero_vE = ex_dict["vib_e_key"]
        else:
            zero_vE = None
        return zero_vE

    def get_reaction_energies(self, zero_vib: bool = True,
                              select_key="min_e_data") -> np.ndarray:
        """
        get the reaction energies (deltaEs) for all asedatabse.

        Returns:
                deltaEs: numpy array
        """
        assert hasattr(self, "_ase_db")
        symbol_key = "chemical_react_symbol"
        select_key = "min_e_data"
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
