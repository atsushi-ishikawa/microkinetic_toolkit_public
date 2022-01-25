#!/usr/bin/env python3


# formal lib
import numpy as np
import pandas as pd
import tempfile
import warnings
import cv2
import os
from typing import List, Iterator, Tuple
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import write as ase_write
from ase.db import connect
from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
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
            w_mes = "set react_e_dict by default state_key and select_key.\n"\
                    "state_key: chem_react_symbol, select_key: min_e_data."
            warnings.warn(w_mes)
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
                elif self.G_option:
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


class ExtraUtilMixin():
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

    def plot_deltaE(self) -> Tuple[Figure, Axes]:
        e_iter = self._gene_reaction_e()
        sym_iter = self._gene_react_symbols()
        fig, axes = plt.subplots()
        y_li = []
        tick_label = []
        for sym, e in enumerate(zip(sym_iter, e_iter)):
            y_li.append(e)
            tick_label.append(sym)
        x_ar = np.arange(len(y_li))
        axes.bar(x_ar, y_li, width=0.5,
                 tick_label=tick_label)
        return fig, axes

    def save_fig_deltaE(self, opng: str):
        fig, axes = self.plot_deltaE()
        fig.savefig(opng)

    @classmethod
    def from_df(cls, df: DataFrame):
        """
        Read elementary reactions from CSV.

        Args:
                csv_file: CSV file with elementary reactions
        Returns:
                Reactions
        """
        reaction_list = []
        for i, row in df.iterrows():
            ddict = row.to_dict()
            reaction = ReactionForDB.from_dict(ddict)
            reaction_list.append(reaction)
        return cls(reaction_list)


def resize_img_to_hmin(im_list):
    h_min = min([im.shape[0] for im in im_list])
    resized_imgs = []
    for im_ins in im_list:
        ini_h, ini_w = im_ins.shape[:2]
        w = round(h_min/ini_h * ini_w)
        new_img = cv2.resize(im_ins,
                             dsize=(w, h_min))
        resized_imgs.append(
                        new_img)
    return resized_imgs


class SnapshotAtomsMixin():
    raise NotImplementedError("")

    def _gene_reactants_fig(self, angles: str = None):
        for reaction in self.reaction_list:
            _ = self._get_atoms_reacts_from_asedb(
                                                     reaction)
        raise NotImplementedError("")

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
                row = rows_list[0]
                atoms = row.toatoms()
                atoms_list.append(atoms)
            elif len(rows_list) == 0:
                warn_mes = "{} doesn't exist in ase database.".format(
                                                              react_sp)
                warnings.warn(
                            warn_mes)
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
                row = rows_list[0]
                atoms = row.toatoms()
                atoms_list.append(atoms)
            elif len(rows_list) == 0:
                warn_mes = "{} doesn't exist in ase database".format(
                                                                react_sp)
                warnings.warn(
                            warn_mes)
            else:
                raise NotImplementedError("")
            return atoms_list

    def _gene_atoms_fig_from_atmslist(self, atoms_list: List[Atoms],
                                      out_dpath: str, reactions_id: int,
                                      angles: str = None):
        "reaction_{id}/{chem_sym}_img.png"
        # snapshotter = SnapShotterAtoms(out_dpath, base_png)
        raise NotImplementedError("")


class SnapShotterAtoms(object):
    def __init__(self, dpath: str,
                 base_png: str,
                 angles=["0z", "10z,-80x", "-90x"]):
        assert os.path.isdir(dpath)
        self.dpath = dpath
        self.base_png = base_png
        self.angles = angles

    def _atoms_to_png(self, atoms: Atoms, opng: str,
                      enlarge_r: tuple = (1, 1, 1),
                      rot_angle: str = "10z,-80x"):
        """
        it writes internal ase atoms object to image file in png format.

        Args:
            opng (str): file name with png extension as a output.
            enlarge_r (tuple): enlarge_r with 3 int objects.
                               default is (1, 1, 1).
                               if you enter (2, 2, 1),
                               2X2X1 supercell are generated as a image.
            rot_angle (str): rot angle designates rotation angle.
                             default "10z,-80x"
        """
        assert opng.endswith(".png")
        ase_write(opng, atoms * enlarge_r,
                  rotation=rot_angle)

    def _atoms_to_cv2img(self, atoms: Atoms,
                         enlarge_r: tuple = (1, 1, 1),
                         rot_angle: str = "10z,-80x") -> np.ndarray:
        fd, path = tempfile.mkstemp(suffix=".png")
        self.atoms_to_png(path, enlarge_r=enlarge_r,
                          rot_angle=rot_angle)
        cv2img = cv2.imread(path)
        os.remove(path)
        return cv2img

    def __call__(self, atoms: Atoms, png_kwds: dict):
        assert isinstance(atoms, Atoms)
        png_fbase = self.base_png.format(
                                    png_kwds)
        png = os.path.join(self.dpath,
                           png_fbase)
        if self.angles is None:
            self._atoms_to_png(png)
        elif len(self.angles) == 1:
            self._atoms_to_png(png,
                               rot_angle=self.angles[0])
        elif len(self.angles) > 1:
            imgs = []
            for rot_angle in self.angles:
                im = self._atoms_to_cv2img(
                                rot_angle=rot_angle)
                imgs.append(im)
            tmp_imgs = resize_img_to_hmin(imgs)
            new_img = cv2.hconcat(tmp_imgs)
            cv2.imwrite(png, new_img)
        else:
            raise AssertionError("")


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
