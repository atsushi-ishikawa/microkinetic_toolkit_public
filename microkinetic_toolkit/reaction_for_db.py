#!/usr/bin/env python3


from .reaction import Reaction
from chem_eq import Reaction as ReactionH


class ReacionForDB(Reaction, ReactionH):
    def get_delta_e(self, e_dict: dict) -> float:
        """set_react_e.

        Args:
            e_dict (dict): e_dict
        """
        lh_E = self.sympy_reactants.subs(
                                     e_dict)
        rh_E = self.sympy_products.subs(
                                     e_dict)
        react_e = rh_E - lh_E
        if not isinstance(react_e, SYMPY_FLOAT):
            emes = "react_e: {}\nit still include Symbol.".format(
                                                                react_e)
            raise AssertionError(emes)
        return float(react_e)

    def _set_sympy_products(self):
        sympy_products = self._get_sympy_hside(
                                        self.products)
        self._sympy_products = sympy_products

    def _set_sympy_reactants(self):
        sympy_reactants = self._get_sympy_hside(
                                        self.reactants)
        self._sympy_reactants = sympy_reactants
    
    def _get_sympy_hside(self,
                         list_nsp: List[N_SPECIE_TYPE]) -> Expr:
        hside_expr = 0.0
        for num, sym in list_nsp:
            if " *" in sym:
                sym = sym.replace("*", "{a}")
            spiecie = Symbol(sym)
            term = num * spiecie
            hside_expr = hside_expr + term
        return hside_expr
