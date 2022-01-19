from ase import Atoms

from ase.db import connect
from ase.optimize import BFGS, FIRE
from ase.io import read, write
from ase.visualize import view

import os
import sys

argvs = sys.argv
dbfile = argvs[1]
entry = argvs[2]

db = connect(dbfile)

try:
    # find molecule by number
    entry = int(entry)
    mol = db.get_atoms(num=entry)
except:
    # find molecule by name
    entry = entry
    mol = db.get_atoms(name=entry)

view(mol)
