#!/usr/bin/env python3


import re
from typing import List, Union, Tuple
from sympy import Symbol
from sympy import Expr
from .reaction import Reaction

INT_NONE_TYPE = Union[int, None]
F_NONE_TYPE = Union[float, None]
N_SPECIE_TYPE = Tuple[int, str]
H_SIDE_TYPE = List[N_SPECIE_TYPE]
SPECIE_TYPE = Union[str, Symbol]


class ReactionForDB(Reaction):
    def __init__(self, reaction_str, reaction_id=None):
        super().__init__(reaction_str, reaction_id=None)
        self._parse_reaction_str_for_sympy()
        self._set_sympy_products()
        self._set_sympy_reactants()

    def _parse_reaction_str_for_sympy(self):
        r_tuples, p_tuples = self._get_reactans_products_for_sympy()
        self.reactants_for_sympy = r_tuples
        self.products_for_sympy = p_tuples

    def _get_reactans_products_for_sympy(self) -> Tuple:
        products = []
        for num, sp, _ in self.products:
            products.append((num, sp))
        reactants = []
        for num, sp, _ in self.reactants:
            reactants.append((num, sp))
        return reactants, products

    def _hside_to_nspli_for_sympy(self, hside_str: str) -> List[N_SPECIE_TYPE]:
        terms = hside_str.split("+")
        list_nspiecie = [self._term_to_num_spiecie_for_sympy(term)
                         for term in terms]
        return list_nspiecie

    def _term_to_num_spiecie_for_sympy(self, term: str) -> N_SPECIE_TYPE:
        term = term.strip()
        reins = re.compile("^[0-9]+")
        re_match = reins.match(term)
        if re_match is None:
            sp = term.strip()
            return (1, sp)
        else:
            num = int(re_match[0])
            sp = reins.sub("", term)
            sp = sp.strip()
            return (num, sp)

    def get_delta_e(self, e_dict: dict) -> float:
        """set_react_e.

        Args:
            e_dict (dict): e_dict
        """
        lh_E = self._sympy_reactants.subs(
            e_dict)
        rh_E = self._sympy_products.subs(
            e_dict)
        try:
            react_e = float(rh_E - lh_E)
        except TypeError:
            import ipdb; ipdb.set_trace()
        return react_e

    def _set_sympy_products(self):
        sympy_products = self._get_sympy_hside(
            self.products_for_sympy)
        self._sympy_products = sympy_products

    def _set_sympy_reactants(self):
        sympy_reactants = self._get_sympy_hside(
            self.reactants_for_sympy)
        self._sympy_reactants = sympy_reactants

    def _get_sympy_hside(self,
                         list_nsp: List[N_SPECIE_TYPE]) -> Expr:
        hside_expr = 0.0
        for num, sym in list_nsp:
            spiecie = Symbol(sym)
            term = num * spiecie
            hside_expr = hside_expr + term
        return hside_expr

    def get_preexponential_test(self, T=300.0, sden=1.0e-5):
        """
        test functions
        """
        import warnings
        warnings.warn("preexponential factor is substituted by test values")
        return 1.0

    def get_atoms_from_moldb(self, spe: str, site: str):
        assert hasattr(self, "_moldb")
        with connect('mols.db') as db:
            query_dict = {}
            raise NotImplementedError("")
            for row in db.connect(**query_dict):
                atoms = row.toatoms()


    def set_moldb(self, moldb: str):
        self._moldb = moldb

    def get_preexponential(self, T=300.0, sden=1.0e-5):
        """
        Calculate pre-exponential factor.
        Note on revserse reaction:
        Pre-exponential factor for reverse reaction is generally not needed
        since reverse pre-exponential is calculated from forward one and equilibrium constant.

        Args:
            T: temperature [K]
            sden: site density [mol/m^2]
        Returns:
            A: pre-exponential factor (float)
        """
        mass_sum = 0
        mass_prod = 1
        rad_all = 0
        d_ave = 0
        rxntype = []

        # calculate mean diameter
        for mol in self.reactants:
            coef, spe, site = mol
            nmol = len(self.reactants)
            # need modification for asedb(reported by oda)
            spe_atom = self.get_atoms_from_adsorbate(spe, site)
            if site == "surf":
                # bare surface
                import ipdb; ipdb.set_trace()
                rxntype.append("surf")
            else:
                #vol = methane.data[mol]['molecular_volume']
                vol = 100.0  # TODO: molecule volume in Angstrom**3
                # molecular radius (in Angstrom) calculated from volume
                rad = 3.0/(4.0*np.pi)*np.cbrt(vol)
                rad *= 10e-10  # Angstrom --> m
                # rad *= 0.182  # do empirical correction based on He, to match the vdW radius

                if nmol == 1 and coef == 2:  # reacion among same species
                    rad_all = 2*rad
                else:
                    rad_all += rad

                # sigma = np.pi*rad_all ** 2  # sigma = pi*(rA + rB)^2
                # rad *= 0.182   # do empirical correction based on He, to match the vdW radius

                d_ave += 2*rad/nmol  # mean diameter

                if spe_atom is None:
                    pass  # bare surface
                else:
                    mass = sum(spe_atom.get_masses())

                mass_sum += mass
                mass_prod *= mass

                if site == "gas":
                    rxntype.append(site)
                else:
                    rxntype.append("surf")

        # end loop over molecule

        if all(rxn == "gas" for rxn in rxntype):
            #
            # gas reaction --- collision theory expression
            #
            red_mass = mass_prod / mass_sum
            red_mass *= amu
            # fac_for  = sigma * np.sqrt(8.0*np.pi*kbolt*T / red_mass) * Nav
            # Eq.3.21 in CHEMKIN Theory manual
            A = Nav*d_ave**2*np.sqrt(8.0*np.pi*kB*T/red_mass)
            # A = 0.5*A  # structural factor
            # A = A*1.0e6  # [m^3] --> [cm^3]
            #
            # sqrt(kbolt*T/mass) [kg*m^2*s^-2*K^-1 * K * kg^-1]^1/2 = [m*s^-1]
            # A = [m^2] * [m*s^-1] * Nav = [m^3*s^-1]*[mol^-1] = [mol^-1*m^3*s^-1]
            # finally A in [mol^-1*cm^3*s^-1]
            #
        elif all(rxn == "surf" for rxn in rxntype):
            #
            # Langmuir-Hinshelwood reaction or desorption --- apply Eyring's transition state theory
            #
            A = kB*T/hplanck/sden  # [s^-1]
            # A = kbolt*T/hplanck # [s^-1]
        else:
            #
            # adsorption --- Hertz-Knudsen or Chemkin
            #
            stick = 0.5
            red_mass = mass_prod/mass_sum
            red_mass *= amu  # [kg]
            #
            # --- Hertz-Knudsen (in concentration form) acturally same with chemkin
            #
            # when having site density information
            # [J*K/kg]^1/2 = [kg*m^2*s^-2/kg]^1/2 = [m*s^-1]
            A = np.sqrt(kB*T/(2.0*np.pi*red_mass))
            # multiplying sticking probability and site density: [m*s^-1] --> [m^3/mol/s]
            A = (stick/sden)*A

        return A
