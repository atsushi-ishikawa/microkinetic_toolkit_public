import re
import numpy as np
from tinydb import Query
from ase.units import create_units
from ase.db import connect

# physical constants
units = create_units('2014')
kB = units['_k']        # Boltzman's constant [J/K]
Nav = units['_Nav']      # Avogadro's constant [mole/mol]
amu = units['_amu']      # atomic mass unig (1.66e-27) [kg]
hplanck = units['_hplanck']  # Plank's constant (6.63e-34) [J*s]
R = kB*Nav             # R (gas constant) [J/mol/K]

# unit conversion
eVtokJ = 96.487
eVtoJ = eVtokJ*1.0e3

# database containing informations for adsorbates
MOLECULE_DB = "molecules.json"


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
        self.products = self._from_halfside_to_nsp_list(products_str)

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
        re_site = re.compile("_(atop|fcc|hcp|br|stable)")
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
            site = re_site_search[1]

        spe = re_site.sub("", re_coef.sub("", term))
        spe = spe.strip()

        # add "_surf" for surface species
        if re_site_search is not None:  # has site
            spe += "_surf"

        return (coef, spe, site)

    def to_dict(self):
        """
        convert to dict.

        Returns:
            base_key (dict)
        """
        def extract_attribute_dict(instance, keys, default=None):
            ddict = {}
            for key in keys:
                val = getattr(instance, key, default)
                ddict[key] = val
            return ddict
        return extract_attribute_dict(self, self.base_keys)

    @classmethod
    def from_dict(cls, ddict):
        reaction_str = ddict["_reaction_str"]
        reaction_id = ddict["_reaction_id"]
        return cls(reaction_str, reaction_id)

    # tinyDB
    def to_tinydb(self, tinydb):
        tinydb.insert(self.to_dict())

    def update_tinydb(self, tinydb):
        ddict = self.to_dict()
        query = Query()
        q_ins = getattr(query, "_reaction_id")
        tinydb.update(ddict, q_ins == ddict["_reaction_id"])

    @property
    def base_keys(self):
        base_keys = ["_reaction_str", "_reaction_id"]
        return base_keys

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

    def get_atoms_from_adsorbate(self, specie, site, surface=None):
        """
        Returns Atoms(H) from H_surf.

        Args:
            specie:
            site:
            surface:
        Returns:
            Atoms
        """
        # MOLECULE_DB is global variables.
        db = connect(MOLECULE_DB)

        if site == "gas":
            if specie == "surf":
                return surface  # already atoms
            else:
                return db.get_atoms(name=specie)
        else:
            spe = specie.split("_")[0]
            return db.get_atoms(name=spe)

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
            db = connect(ase_db)
            formula = atom.get_chemical_formula()
            try:
                energy = db.get(name=formula).data.energy
            except:
                print("{0:} is not found in {1:}.".format(formula, db))
                energy = 0.0

            return energy

        def adsorb_on_surface(ads=None, surface=None, height=2.0):
            """
            Adsorb molecule on surface.

            Args:
                ads:
                surface:
                height:
            Returns:
                atoms
            """
            from ase.build import add_adsorbate
            from ase.visualize.plot import plot_atoms
            from ase.io import write
            from .preparation.rot_control import make_atoms_with_standard_alignment
            import matplotlib.pyplot as plt

            #surf_copy = copy.deepcopy(surface)  # deepcopy does not work with Cython
            surf_copy = surface.copy()
            surf_copy.calc = surface.calc

            # adjust height
            shift = min(ads.positions[:, 2])
            height -= shift

            # TODO: rotate atom
            #ads = make_atoms_with_standard_alignment(ads)

            add_adsorbate(surf_copy, ads, offset=(0, 0), position=(0, 0), height=height)
            surf_copy.pbc = True

            poscar = True
            snapshot = False

            if snapshot:
                fig, ax = plt.subplots()
                plot_atoms(surf_copy)
                ax.set_axis_off()
                fig.savefig(surf_copy.get_chemical_formula() + ".png")
                plt.close()
            if poscar:
                filename = surf_copy.get_chemical_formula() + "_POSCAR"
                write(filename, surf_copy)

            return surf_copy

        def optimize_geometry(atoms):
            import ase.optimize

            #atoms_copy = copy.deepcopy(atoms)  # deepcopy does not work with Cython
            atoms_copy = atoms.copy()
            atoms_copy.calc = atoms.calc

            opt = ase.optimize.bfgs.BFGS(atoms_copy, maxstep=0.1)
            opt.run(steps=30)

            return atoms_copy

        # main
        # set calculator
        if calculator == "vasp":
            # TODO: do vasp setting
            pass
        elif calculator == "emt":
            calc = EMT()
        else:
            raise Exception("use vasp or emt for calculator")

        energy_dict = {"reactants": 0.0, "products": 0.0}
        for side in ["reactants", "products"]:
            sequence = self.reactants if side == "reactants" else self.products
            for i in sequence:
                _, mol, site = i
                atoms = self.get_atoms_from_adsorbate(mol, site, surface=surface)

                if site != "gas":
                    atoms = adsorb_on_surface(ads=atoms, surface=surface, height=1.5)

                # look up ase database
                if ase_db is not None:
                    energy = search_energy_from_ase_db(atoms, ase_db)
                else:
                    atoms.set_calculator(calc)
                    atoms = optimize_geometry(atoms)
                    energy = atoms.get_potential_energy()

                energy_dict[side] += energy

        deltaE = energy_dict["products"] - energy_dict["reactants"]

        return deltaE

    def get_entropy_difference(self):
        """
        Calculate entropy difference (deltaS, in eV/K) along the reaction.
        Moleclar entropy should be prepared beforehand.

        Returns:
            deltaS (float)
        """
        db = connect(MOLECULE_DB)

        entropies = {"reactants": 0.0, "products": 0.0}
        for side in ["reactants", "products"]:
            sequence = getattr(self, side)
            for mol in sequence:
                _, spe, site = mol

                if site != "gas" or spe == "surf":
                    # surface species
                    entropy = 0.0
                else:
                    try:
                        atoms = db.get(name=spe)
                        entropy = atoms.key_value_pairs["molecular_entropy"]  # [eV/K]
                    except:
                        entropy = 1.0e-3  # rough entropy estimation ... 1 meV/K

                entropies[side] += entropy

        deltaS = entropies["products"] - entropies["reactants"]
        return deltaS

    def get_deltaN(self):
        """
        Get difference in the gas-phase coefficient of reactants and products.
        The value has contribution on deltaG as P*V*deltaN.

        Returns:
            deltaN (float)
        """
        gas_coefs = {"reactants": 0.0, "products": 0.0}
        for side in ["reactants", "products"]:
            sequence = getattr(self, side)
            for mol in sequence:
                coef, spe, site = mol
                if site == "gas" and spe != "surf":
                    # gas-phase species
                    gas_coefs[side] += coef

        deltaN = gas_coefs["products"] - gas_coefs["reactants"]
        return deltaN

    def get_deltaV(self):
        """
        Get difference in the gas-phase volume of reactants and products.
        The value has contribution on deltaG as P*deltaV*deltaN.

        Returns:
            deltaV (float) [m^3/mol]
        """
        db = connect(MOLECULE_DB)

        V = {"reactants": 0.0, "products": 0.0}
        for side in ["reactants", "products"]:
            sequence = getattr(self, side)
            for mol in sequence:
                _, spe, site = mol

                if site != "gas" or spe == "surf":
                    # surface species
                    volume = 0.0
                else:
                    try:
                        atoms   = db.get(name=spe)
                        volume  = atoms.key_value_pairs["molecular_volume"]  # Angstrom^3
                        volume *= 1.0e-30*Nav  # Angstrom^3 --> m^3/mol
                    except:
                        volume  = 100.0  # rough volume estimation ... 100 Angstrom^3
                        volume *= 1.0e-30*Nav  # Angstrom^3 --> m^3/mol

                V[side] += volume

        deltaV = V["products"] - V["reactants"]
        return deltaV

    def get_rate_constant(self, deltaE=None, T=300.0, bep_param: dict = None, sden=1.0e-5):
        """
        Calculate rate constant from reaction energy (deltaE).

        Args:
            deltaE: reaction energy [eV]
            T: temperature [K]
            bep_param: BEP alpha and beta (for eV unit)
            sden: site density [mol/m^2]
        Returns:
            rate constant
        """
        # calculate Ea from deltaE
        Ea = bep_param["alpha"]*deltaE + bep_param["beta"]
        Ea *= eVtoJ

        A = self.get_preexponential(T=T, sden=sden)

        exp = np.exp(-Ea/R/T)
        rateconst = A*exp

        return rateconst

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

            if spe == "surf":
                # bare surface
                rxntype.append("surf")
            else:
                #vol = methane.data[mol]['molecular_volume']
                vol = 100.0  # TODO: molecule volume in Angstrom**3
                # molecular radius (in Angstrom) calculated from volume
                rad = 3.0/(4.0*np.pi)*np.cbrt(vol)
                rad *= 1.0e-10  # Angstrom --> m
                # rad *= 0.182  # do empirical correction based on He, to match the vdW radius

                if nmol == 1 and coef == 2:  # reacion among same species
                    rad_all = 2*rad
                else:
                    rad_all += rad

                # sigma = np.pi*rad_all ** 2  # sigma = pi*(rA + rB)^2
                # rad *= 0.182   # do empirical correction based on He, to match the vdW radius

                d_ave += 2*rad/nmol  # mean diameter

                if spe_atom is None:
                    # spe_atom in only molecules.
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
            with open("gas_react1.txt", "a") as add:
                print(self.reactants, file=add)
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
            with open("langmuir_or_desorpton2.txt", "a") as add:
                print(self.reactants, file=add)
            #
            # Langmuir-Hinshelwood reaction or desorption --- apply Eyring's transition state theory
            #
            A = kB*T/hplanck/sden  # [s^-1]
        else:
            with open("adsorption_react3.txt", "a") as add:
                print(self.reactants, file=add)
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
