import numpy as np
from .trn import *
from . import cvt

from chimerax.core.commands import run




class PointSampler:


	def __init__(self, volume, random_seed=None):
		self.volume = volume
		self.x = []
		self.y = []
		self.z = []
		self.random_seed = random_seed
		self.map_th = self.volume.data.matrix().copy()

	def threshold(self, thresh):
		self.map_th[self.map_th < thresh] = 0

	def sample(self, M):
		pass

	def get_global_coords(self):
		if self.volume is None:
			return self.x, self.y, self.z

		global_x = []
		global_y = []
		global_z = []
		for i in range(len(self.x)):
			coord = self.volume.ijk_to_global_xyz([self.z[i],self.y[i],self.x[i]])
			global_x.append(coord[0])
			global_y.append(coord[1])
			global_z.append(coord[2])
		return global_x, global_y, global_z

	def get_global_centered_coords(self):
		global_x, global_y, global_z = self.get_global_coords()
		centroid = self.get_centroid()
		for i in range(len(self.x)):
			global_x[i] -= centroid[0]
			global_y[i] -= centroid[1]
			global_z[i] -= centroid[2]
		return global_x, global_y, global_z

	def show_points(self, max_points=None, size=2, color='#ff5a5a'):
		if max_points is None:
			max_points = len(self.x)

		global_x, global_y, global_z = self.get_global_coords()
		for i in range(max_points):
			run(self.volume.session, "shape sphere radius %d center %f,%f,%f color %s"%(size, global_x[i], global_y[i], global_z[i], color))

	def save_points(self, fname):
		global_x, global_y, global_z = self.get_global_coords()
		f = open(fname, 'w')
		f.write(json.dumps({'x':global_x, 'y':global_y, 'z':global_z}))
		f.close()

	def load_points(self, fname):
		f = open(fname, 'r')
		j = json.loads(f.read())
		f.close()
		self.volume = None
		self.x = j['x']
		self.y = j['y']
		self.z = j['z']

	def get_centroid(self):
		centroid = np.array([0., 0., 0.])
		global_x, global_y, global_z = self.get_global_coords()
		for i in range(len(global_x)):
			centroid[0] += global_x[i] / len(global_x)
			centroid[1] += global_y[i] / len(global_x)
			centroid[2] += global_z[i] / len(global_x)
		return centroid





class TRNPointSampler(PointSampler):


	def __init__(self, volume, random_seed=None, l0_factor=0.005, lf=0.5, tf_factor=8, e0=0.3, ef=0.05):
		super().__init__(volume, random_seed)
		self.l0_factor = l0_factor
		self.lf = lf
		self.tf_factor = tf_factor
		self.e0 = e0
		self.ef = ef

	def sample(self, M):
		rm0,arr_flat,arr_idx,xyz,coords_1d = trn_rm0(self.map_th, M, random_seed=self.random_seed)
		l0 = self.l0_factor * M 
		lf = self.lf
		tf = M * self.tf_factor
		e0 = self.e0
		ef = self.ef

		rms,rs,ts_save = trn_iterate(rm0,arr_flat,arr_idx,xyz,n_save=10,e0=e0,ef=ef,l0=l0,lf=lf,tf=tf,do_log=True,log_n=10)

		N_cube = max(self.map_th.shape[0],self.map_th.shape[1],self.map_th.shape[2])
		N_cube += N_cube%2
		
		for p in rms[10]:
			self.x.append(p[0] + N_cube//2)
			self.y.append(p[1] + N_cube//2)
			self.z.append(p[2] + N_cube//2)





class CVTPointSampler(PointSampler):


	def __init__(self, volume, random_seed=None, max_iter=10):
		super().__init__(volume, random_seed)
		self.max_iter = max_iter

	def sample(self, M):
		cvt.session = self.volume.session

		robs, map_edited = cvt.get_init(self.map_th, M, self.random_seed)
		# self.volume.session.logger.info(str(list(robs)))

		robs = cvt.iterate(map_edited, robs)

		N = map_edited.shape[0]
		for i in range(len(robs)):
			self.x.append(robs[i][0] * N)
			self.y.append(robs[i][1] * N)
			self.z.append(robs[i][2] * N)
			

	

	