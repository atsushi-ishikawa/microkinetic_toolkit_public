import numpy as np

def calculate_volume_and_entropy():
	"""
	Calculate volume and entropy for all the molecules in the current reaction set.

	Returns:
	"""
	from ase import Atoms, Atom
	from ase.calculators.gaussian import Gaussian
	from ase.db import connect
	from ase.optimize import BFGS
	from ase.vibrations import Infrared

	calculator = "gau"
	calculator = calculator.lower()

	oldjson = "tmp.json"
	newjson = "tmp_new.json"

	db1 = connect(oldjson)
	db2 = connect(newjson)

	add_volume = True
	add_entropy = True

	# Gaussian
	if "gau" in calculator:
		method = "b3lyp"
		basis = "6-31G*"
	else:
		raise RuntimeError('molecular volume only implemented in Gaussian')

	# add volume and entropy information for all the molecules in database
	Nmol = len(db1)
	for imol in range(Nmol):  # add 10 for safety
		try:
			id = imol + 1
			mol = db1.get(id=id).name
			tmp = db1.get_atoms(id=id)
			print("adding volume and entropy", mol)

			# do geometry optimization first
			tmp.calc = Gaussian(method=method, basis=basis)
			BFGS(tmp).run(fmax=0.05)

			# volume and molecular total entropy calculation
			if add_volume and not add_entropy:
				tmp.calc = Gaussian(method=method, basis=basis, volume='tight')
				tmp.get_potential_energy()
				vol = tmp.get_molecular_volume()
			if not add_volume and add_entropy:
				tmp.calc = Gaussian(method=method, basis=basis,
					force=None, freq='noraman')
				tmp.get_potential_energy()
				entropy = tmp.get_molecular_entropy()
			if add_volume and add_entropy:
				tmp.calc = Gaussian(method=method, basis=basis, volume='tight',
					force=None, freq='noraman')
				tmp.get_potential_energy()
				vol = tmp.get_molecular_volume()
				entropy = tmp.get_molecular_entropy()

			# look for magmom
			try:
				magmom = tmp.get_magnetic_moments()
			except:
				magmom = tmp.get_initial_magnetic_moments()

			tmp.set_initial_magnetic_moments(magmom)

			# now write to database
			if add_volume and add_entropy:
				db2.write(tmp, key_value_pairs={'name': mol,
												'molecular_volume': vol,
												'molecular_entropy': entropy})
			elif add_volume and not add_entropy:
				db2.write(tmp, key_value_pairs={'name': mol,
												'molecular_volume': vol})
			elif not add_volume and add_entropy:
				db2.write(tmp, key_value_pairs={'name': mol,
												'molecular_entropy': entropy})
		except:
			print("has nothing --- go to next")

	return None

def make_atoms_standard_alignment(atoms):
	"""
	Make atoms to standard alignment.
	Taken from oda-san's make_atoms_with_stdpos of lib/rot_control.py.

	Args:
		atoms:
	Returns:
		atoms:
	"""
	from numpy import linalg

	w_cov = np.cov(atoms.positions.T, bias=True, aweights=atoms.get_masses())
	w_, v_ = linalg.eig(w_cov)
	ids = np.argsort(w_)[::-1]

	# the max variance of atoms is along z-axis.
	v_ = v_[ids]
	T_mat = linalg.inv(v_)
	cp_atoms = atoms.copy()
	new_pos = np.dot(cp_atoms.positions, T_mat)
	cp_atoms.positions = new_pos
	return cp_atoms
