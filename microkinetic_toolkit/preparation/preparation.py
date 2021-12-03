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

def make_atoms_with_standard_alignment(atoms, direction=None):
	"""
	Make atoms to standard alignment.
	Copied from oda-sans make_atoms_with_stdpos in libs/rot_control.py

	Args:
		atoms:
	Returns:
		atoms:
	"""

	w_cov = np.cov(atoms.positions.T, bias=True, aweights=atoms.get_masses())
	w_, v_ = np.linalg.eig(w_cov)
	if direction == "side":
		ids = np.argsort(w_)[::-1]  # ::-1 means reverse
	else:
		ids = np.argsort(w_)

	# the max variance of atoms is along z-axis.
	v_ = v_[ids]
	T_mat = np.linalg.inv(v_)
	atoms_copy = atoms.copy()
	new_pos = np.dot(atoms_copy.positions, T_mat)
	atoms_copy.positions = new_pos

	return atoms_copy

def rotate_atoms_accoridngto_rotvec(atoms, rotvec):
	atoms_copy = atoms.copy()
	rad = np.linalg.norm(rotvec)
	deg_angle = np.rad2deg(rad)
	if np.isclose(rad, 0):
		pass
	else:
		rotaxis = rotvec / rad
		com = atom_copy.get_center_of_mass()
		atoms_copy.rotate(deg_angle, rotaxis, center=com)
	return atoms_copy

def get_qternion_from_v1tov2(v1, v2):
	rotaxis = np.cross(v1, v2)
	rotaxis = rot_axis / np.linalg.norm(rotaxis)
	rad = calc_angle(v1, v2)
	rotvec = rad*rotaxis
	qternion = quaternion.from_rotation_vector(rotvec)
	return qternion

def get_rotvec_from_v1tov2(v1, v2):
	qternion = get_qternion_from_v1tov2(v1, v2)
	rotvec = quaternion.as_rotation_vector(qternion)
	return rotvec

def get_rotvec_into_minus_z(atoms, atom_id=0, rotaxis=None):
	target_atom_pos = atoms.positions[atom_id]
	com = atoms.get_center_of_mass()
	ref_vec = target_atom_pos - com
	minus_z = np.array([0, 0, -1])
	if rotaxis is None:
		rotaxis = np.cross(ref_vec, minus_z)
		norm = np.linalg.norm(rotaxis)
		if np.isclose(norm, 0):
			rotaxis = np.array([1.0, 0.0, 0.0])
		else:
			rotaxis = rotaxis / np.linalg.norm(rotaxis)
		rad = calc_angle(ref_vec, mi_z)
		rotvec = rad * rotaxis
	else:
		orthogonal_vec = get_orthogonal_vec(ref_vec, rotaxis)
		n = np.linalg.norm(orthogonal_vec)
		if np.isclose(n, 0):
			rotvec = np.array([0, 0, 0])
		else:
			norm_vec = get_norm_vec(orthogonal_vec)
			rotvec = get_rotvec_from_v1tov2(norm_vec, mi_z)
	return rotvec

def rotate_atoms_into_minus_z(atoms, atoms_id=0, rotaxis=None):
	rotvec = get_rotvec_into_minus_z(atoms, atoms_id=atoms_id, rotaxis=rotaxis)
	new_atoms = rotate_atoms_accoridngto_rotvec(atoms, rotvec=rotvec)
	return new_atoms

def prepare(species):
	from ase import Atoms
	from ase.build import molecule
	from ase.visualize import view

	for mol in species:
		print(mol)
		atoms = Atoms(molecule(mol))
		print(atoms.positions)
		#view(atoms)
		rotated = make_atoms_with_standard_alignment(atoms, direction="side")
		print(rotated.positions)
		view(rotated)
	quit()
	return None
