#!/usr/bin/env python3

import numpy as np
import argparse
from ase.db import connect

def calculate_volume_and_entropy(Atoms=None, names=None, json=None, calculator="emt", vibration_freq=None):
    """
    Calculate volume and entropy for all the molecules in the current reaction set.
    Assuming that,
      1) geometry optimization is already done
      2) vibrational analysis is already done thus frequencies are given (otherwise vib. entropy is neglected)
    Args:
        Atoms:
        names:
        json:
        calculator:
        vibration_freq:
    Returns:
        None
    """
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.optimize import QuasiNewton
    from ase.vibrations import Vibrations
    from ase.thermochemistry import IdealGasThermo
    from ase.visualize import view

    calculator = calculator.lower()
    db = connect(json)

    if not isinstance(names, list):
        names = [names]

    atoms = Atoms
    atoms.calc = EMT()


    if vibration_freq is None:
        vib_energies = 3 * len(atoms) * [1.0e1]
    else:
        vib_energies = vibration_freq

    linear_molecules = ["CO2", "C2H2", "C2H", "HCN2", "HCNO", "NO", "ON"]
    if len(atoms) == 1:
        geometry = "monatomic"
    elif len(atoms) == 2 or str(atoms.symbols) in linear_molecules:
        geometry = "linear"
    else:
        geometry = "nonlinear"

    print(" ---{}--- {}".format(str(atoms.symbols), geometry))

    potentialenergy = atoms.get_potential_energy()
    thermo = IdealGasThermo(atoms=atoms, vib_energies=vib_energies, potentialenergy=potentialenergy,
                            geometry=geometry, symmetrynumber=1, spin=0)
    entropy = thermo.get_entropy(temperature=298.15, pressure=101325.0)
    #
    # volume
    #
    if calculator == "gaussian":
        method = "b3lyp"
        basis  = "6-31G*"
        atoms.calc = Gaussian(method=method, basis=basis)
        atoms.get_potential_energy()
        volume = atoms.get_molecular_volume()
    else:
        vdw_radii = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52}  # Bondi [Angstrom]
        volume = 0  # molecular volume [Angstrom^3]
        for atom in atoms:
            radii = vdw_radii[atom.symbol]
            volume += 4.0 / 3.0 * np.pi * radii ** 3
    #
    # add to database
    #
    for name in names:
        key_value_pairs = {"molecular_entropy": entropy, "name": name, "molecular_volume": volume}
        db.write(atoms, key_value_pairs=key_value_pairs)

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", default=None, help="name of json file to add data")
    parser.add_argument("--database", default=None, help="ASE/json database containing Atoms object")
    parser.add_argument("--names", default=None, help="list of name: can set multiple name to the same Atoms")
    parser.add_argument("--num", default="all", help="number of records (or id) in database, to add to \n"
                                                     "new json file. If add all the records, set \"all\"")
    args = parser.parse_args()
    json = args.json
    database = args.database
    names = args.names
    num = args.num

    db = connect(database)

    if num == "all":
        range = range(0, len(db))
    else:
        range = range(num, num+1)

    for id in range:
        try:
            atoms = db.get_atoms(id=id)
            if names is None:
                names_ = str(atoms.symbols)
            else:
                names_ = names
            calculate_volume_and_entropy(Atoms=atoms, names=names_, json=json)
        except:
            continue
