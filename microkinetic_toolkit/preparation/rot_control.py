import numpy as np

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

def get_rotvec_from_v1tov2(v1, v2):
	# waiting for oda-san
	return None

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

def make_atoms_with_standard_alignment(atoms, direction=None):
	"""
	Make atoms to standard alignment.
	Copied from oda-sans make_atoms_with_stdpos in libs/rot_control.py

	Args:
		atoms:
		direction:
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
