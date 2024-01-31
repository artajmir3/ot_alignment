import math
import numpy as np

from scipy.spatial.transform import Rotation

from chimerax.core.commands import run
from chimerax.geometry import Place

class Term:
	def __init__(self, coef=1, exps=[0,0,0,0]):
		self.coef = coef
		self.exps = exps
		while len(self.exps) < 8:
			self.exps.append(0)
	
	def derivative(self, index):
		new_exps = [0,0,0,0,0,0,0,0]
		for i in range(8):
			if i != index:
				new_exps[i] = self.exps[i]
			else:
				new_exps[i] = self.exps[i] - 1
		return Term(coef=self.coef*self.exps[index], exps=new_exps)
	
	def evaluate(self, vals):
		s = self.coef
		if s == 0:
			return 0
		for i in range(min(len(vals), len(self.exps))):
			if self.exps[i] == 0:
				continue
			elif self.exps[i] == 1:
				s *= vals[i]
			elif self.exps[i] == 2:
				s *= vals[i] * vals[i]
			elif self.exps[i] == 3:
				s *= vals[i] * vals[i] * vals[i]
			else:
				s *= vals[i] ** self.exps[i]
		return s
	
	def bunch_evaluate(self, bunch_vals):
		s = self.coef
		for i in range(min(len(bunch_vals), len(self.exps))):
			for j in range(self.exps[i]):
				s = np.multiply(s, bunch_vals[i])
		return s
	
	def __mul__(self, other):
		coef = self.coef * other.coef
		exps = []
		for i in range(8):
			exps.append(self.exps[i] + other.exps[i])
		return Term(coef=coef, exps=exps)
	
	def __neg__(self):
		return Term(coef=-self.coef, exps=self.exps)
	
	def __str__(self):
		s = ""
		for i in range(8):
			if self.exps[i] != 0:
				s += "q" + str(i) + "^" + str(self.exps[i])
		return str(self.coef) + "*" + s
	
	def is_similar(self, other):
		for i in range(8):
			if self.exps[i] != other.exps[i]:
				return False
		return True
	
	def simplify(self, vals):
		new_exps = [0,0,0,0,0,0,0,0]
		coef = self.coef
		for i in range(8):
			if vals[i] != 'x':
				coef *= vals[i]**self.exps[i]
			else:
				new_exps[i] = self.exps[i]
		return Term(coef=self.coef*self.exps[index], exps=new_exps)
		
	
class Polynomial:
	def __init__(self, terms=None):
		if terms is None:
			terms = []
		self.terms = {}
		for term in terms:
			self.terms[tuple(term.exps)] = term
		
	def add_term(self, term):
		if term.coef == 0:
			return
		if tuple(term.exps) not in self.terms:
			self.terms[tuple(term.exps)] = term
		else:
			self.terms[tuple(term.exps)].coef += term.coef

	def get_terms(self):
		terms = []
		for exps in self.terms:
			terms.append(self.terms[exps])
		return terms
		
	def derivative(self, index):
		res = Polynomial()
		for term in self.get_terms():
			res.add_term(term.derivative(index))
		return res
	
	def simplify(self, vals):
		res = ()
		for term in self.get_terms():
			res.add_term(term.simplify(vals))
		return res
	
	def evaluate(self, vals):
		t = 0
		for term in self.get_terms():
			t += term.evaluate(vals)
		return t
	
	def bunch_evaluate(self, bunch_vals):
		t = 0
		for term in self.get_terms():
			t = np.add(t, term.bunch_evaluate(bunch_vals))
		return t
	
	def __add__(self, other):
		res = Polynomial()
		for term in self.get_terms():
			res.add_term(term)
		for term in other.get_terms():
			res.add_term(term)
		return res
	
	def __sub__(self, other):
		res = Polynomial()
		for term in self.get_terms():
			res.add_term(term)
		for term in other.get_terms():
			res.add_term(-term)
		return res
	
	def __neg__(self):
		res = Polynomial()
		for term in self.get_terms():
			res.add_term(-term)
		return res
	
	def __mul__(self, other):
		res = Polynomial()
		for term1 in self.get_terms():
			for term2 in other.get_terms():
				res.add_term(term1*term2)
		return res
	
	def __str__(self):
		s = ""
		terms = self.get_terms()
		for i in range(len(terms)):
			if i > 0:
				s += " + "
			s += str(terms[i])
		return s

class Quaternion:
	def __init__(self, real_pol, i_pol, j_pol, k_pol):
		self.real_pol = real_pol
		self.i_pol = i_pol
		self.j_pol = j_pol
		self.k_pol = k_pol
		
	def conjugate(self):
		return Quaternion(self.real_pol, -self.i_pol, -self.j_pol, -self.k_pol)
	
	def __mul__(self, other):
		return Quaternion(self.real_pol*other.real_pol - self.i_pol*other.i_pol - self.j_pol*other.j_pol - self.k_pol*other.k_pol,
						  self.real_pol*other.i_pol + self.i_pol*other.real_pol + self.j_pol*other.k_pol - self.k_pol*other.j_pol,
						  self.real_pol*other.j_pol + self.j_pol*other.real_pol + self.k_pol*other.i_pol - self.i_pol*other.k_pol,
						  self.real_pol*other.k_pol + self.k_pol*other.real_pol + self.i_pol*other.j_pol - self.j_pol*other.i_pol)
	
	def __str__(self):
		return str(self.real_pol) + " + (" + str(self.i_pol) + ")i + (" + str(self.j_pol) + ")j + (" + str(self.k_pol) + ")k"


q = Quaternion(Polynomial([Term(coef=1, exps=[1,0,0,0])]), Polynomial([Term(coef=1, exps=[0,1,0,0])]), 
			   Polynomial([Term(coef=1, exps=[0,0,1,0])]), Polynomial([Term(coef=1, exps=[0,0,0,1])]))
qs = q.conjugate()
p = Quaternion(Polynomial([Term(coef=0, exps=[0,0,0,0])]), Polynomial([Term(coef=1, exps=[0,0,0,0])]), 
               Polynomial([Term(coef=1, exps=[0,0,0,0])]), Polynomial([Term(coef=1, exps=[0,0,0,0])]))

b = q*p*qs


def get_quaternion_vals(theta, ax, ay, az):
	"""
	Compute the quaternion representation for a given rotation in angle-axis representation
	params:
		theta: the angle of the rotation in radians
		ax, ay, az: three floats in a way that *ax, ay, az) shows the 3d axis of the rotation

	retrun:
		q: is a list of length 4 that has values of the corresponding quaternion
	"""
	
	n = math.sqrt(ax**2 + ay**2 + az**2)
	return [math.cos(theta/2), math.sin(theta/2)*ax/n, math.sin(theta/2)*ay/n, math.sin(theta/2)*az/n]

def convert_to_poly(vals):
	return Quaternion(Polynomial([Term(coef=vals[0], exps=[0,0,0,0])]), Polynomial([Term(coef=vals[1], exps=[0,0,0,0])]), 
					  Polynomial([Term(coef=vals[2], exps=[0,0,0,0])]), Polynomial([Term(coef=vals[3], exps=[0,0,0,0])]))
	

def perform(x, y, z, vals):
	"""
	Apply a given rotation on a given point cloud and generate a new point cloud
	params:
		x, y, z: three lists with len(x)=len(y)=len(z) in a way that (x[i], y[i], z[i]) is the 3d coordinates of the i-th point
		vals: a list of length 4 that contains the values of the quaternion correponding the the ritation

	return:
		xr, yr, zr: three lists with len(xr)=len(yr)=len(zr) in a way that (xr[i], yr[i], zr[i]) is the 3d coordinates of the i-th point after the rotation
	"""
	
	xr = []
	yr = []
	zr = []
	q = convert_to_poly(vals)
	qs = q.conjugate()
	p = Quaternion(Polynomial([Term(coef=1, exps=[1,0,0,0])]), Polynomial([Term(coef=1, exps=[0,1,0,0])]), 
				   Polynomial([Term(coef=1, exps=[0,0,1,0])]), Polynomial([Term(coef=1, exps=[0,0,0,1])]))
	b = q*p*qs
	
	bunch_vals = [np.zeros(len(x)), np.array(x), np.array(y), np.array(z)]
	xr = list(b.i_pol.bunch_evaluate(bunch_vals))
	yr = list(b.j_pol.bunch_evaluate(bunch_vals))
	zr = list(b.k_pol.bunch_evaluate(bunch_vals))
	
	return xr, yr, zr
	
def change_quat_format(vals):
	return [vals[3], vals[0], vals[1], vals[2]]

def diff_quaternions(q1, q2):
	diff = (convert_to_poly(q1) * convert_to_poly(q2).conjugate()).real_pol.evaluate((0,0,0,0))
	if diff > 1:
		diff = 1
	elif diff < -1:
		diff = -1
	deg = math.acos(diff) * 2 * 180 / math.pi
	if deg > 180:
		return 360 - deg
	else:
		return deg

# def transform_map(volume, alignment):
# 	v = volume.writable_copy(require_copy = True, copy_colors = False)
# 	volume.session.logger.info(str(v.id))
# 	volume.session.logger.info(str(volume.id))
# 	Bbar = alignment['Bbar']
# 	Abar = alignment['Abar']
# 	rot = None
# 	if 'q' not in alignment.keys():
# 		rot = change_quat_format(alignment['r'].as_quat())
# 	else:
# 		rot = alignment['q']
# 	volume.session.logger.info(str(rot))
	
# 	volume.session.logger.info(str(v.ijk_to_global_xyz([0,0,0])))

# 	run(v.session, "move %f,%f,%f models #%d"%(-Bbar[0], -Bbar[1], -Bbar[2], v.id[0]))
# 	volume.session.logger.info(str(v.ijk_to_global_xyz([0,0,0])))

# 	v = v.writable_copy(require_copy = True, copy_colors = False)

# 	run(v.session, "turn %f,%f,%f %f center 0,0,0 models #%d"%(rot[1], rot[2], rot[3], math.acos(rot[0]) * 180 / math.pi *2, v.id[0]))

# 	v = v.writable_copy(require_copy = True, copy_colors = False)

# 	run(v.session, "move %f,%f,%f models #%d"%(Abar[0], Abar[1], Abar[2], v.id[0]))

def transform_map(volume, alignment):
	v = volume.writable_copy(require_copy = True, copy_colors = False)
	Bbar = alignment['Bbar']
	Abar = alignment['Abar']
	rot = None
	if 'q' not in alignment.keys():
		rot = change_quat_format(alignment['r'].as_quat())
	else:
		rot = alignment['q']
	volume.session.logger.info(str(rot))

	x = [1,0,0]
	y = [0,1,0]
	z = [0,0,1]
	xr,yr,zr = perform(x, y, z, rot)
	axes_rot = np.asarray([[xr[0], yr[0], zr[0]], [xr[1], yr[1], zr[1]], [xr[2], yr[2], zr[2]]])
	axes_nor = np.asarray([[1,0,0],[0,1,0],[0,0,1]])
	origin_not = np.asarray([0, 0, 0])
	origin_Bbar = np.asarray([-Bbar[0], -Bbar[1], -Bbar[2]])
	origin_Abar = np.asarray([Abar[0], Abar[1], Abar[2]])

	p_Bbar = Place(axes=axes_nor, origin=origin_Bbar)
	p_Abar = Place(axes=axes_nor, origin=origin_Abar)
	p_rot = Place(axes=axes_rot, origin=origin_not)
	# p = Place(axes=np.asarray([[1,0,0],[0,1,0],[0,0,1]]), origin=np.asarray([10,0,0]))

	# run(v.session, "move %f,%f,%f models #%d"%(-Bbar[0], -Bbar[1], -Bbar[2], v.id[0]))
	v.scene_position = p_Bbar * v.scene_position
	# volume.session.logger.info(str(v.ijk_to_global_xyz([0,0,0])))

	v = v.writable_copy(require_copy = True, copy_colors = False)

	v.scene_position = p_rot * v.scene_position
	# run(v.session, "turn %f,%f,%f %f center 0,0,0 models #%d"%(rot[1], rot[2], rot[3], math.acos(rot[0]) * 180 / math.pi *2, v.id[0]))

	v = v.writable_copy(require_copy = True, copy_colors = False)

	v.scene_position = p_Abar * v.scene_position
	# run(v.session, "move %f,%f,%f models #%d"%(Abar[0], Abar[1], Abar[2], v.id[0]))