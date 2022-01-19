import os
from reaction_tools import read_speciesfile
from reaction_tools import remove_parentheses
from ase.collections import methane

massfile = "mass.txt"  # write
species = read_speciesfile("species.txt")
fmass = open(massfile, "w")

masslist = []

for i in species:
    mol = methane[i]
    mass = sum(mol.get_masses())

    masslist.append(mass)

fmass.write(str(masslist))
fmass.close()

remove_parentheses(massfile)
