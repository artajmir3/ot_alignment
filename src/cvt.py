import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.spatial
import sys
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from scipy.spatial import ConvexHull, Delaunay
import time

from .coords import *


# eps = sys.float_info.epsilon
eps = 0.000001
session = None
time_points_generation = 0
time_in_computation = 0
time_point_computation = 0

def in_box(robots, bounding_box):
	return np.logical_and(np.logical_and(bounding_box[0] <= robots[:, 0],
									 robots[:, 0] <= bounding_box[1]),
					  np.logical_and(bounding_box[2] <= robots[:, 1],
									 robots[:, 1] <= bounding_box[3]),
					  np.logical_and(bounding_box[4] <= robots[:, 2],
									 robots[:, 2] <= bounding_box[5]))


def voronoi(robots, bounding_box):
	i = in_box(robots, bounding_box)
	points_center = robots[i, :]
	points_left = np.copy(points_center)
	points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])
	points_right = np.copy(points_center)
	points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])
	points_down = np.copy(points_center)
	points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])
	points_up = np.copy(points_center)
	points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])
	points_back = np.copy(points_center)
	points_back[:, 2] = bounding_box[4] - (points_back[:, 2] - bounding_box[4])
	points_forth = np.copy(points_center)
	points_forth[:, 2] = bounding_box[5] + (bounding_box[5] - points_forth[:, 2])
	points = np.append(points_center,
					   np.append(np.append(points_left,
										   points_right,
										   axis=0),
								 np.append(np.append(points_down,
													   points_up,
													   axis=0),
										   np.append(points_back,
													   points_forth,
													   axis=0),
										   axis=0),
								 axis=0),
					   axis=0)
	#TODO
	# return points
	# Compute Voronoi
	vor = sp.spatial.Voronoi(points)
	# Filter regions and select corresponding points
	regions = []
	points_to_filter = [] # we'll need to gather points too
	ind = np.arange(points.shape[0])
	ind = np.expand_dims(ind,axis= 1)

	# session.logger.info(str(list(vor.points)))
	# session.logger.info(str(list(vor.regions)))
	# session.logger.info(str(list(vor.point_region)))

	for i,region in enumerate(vor.regions): # enumerate the regions
		if not region: # nicer to skip the empty region altogether
			continue

		flag = True
		# session.logger.info(str(bounding_box))
		# session.logger.info('New region %d'%(i,))
		tot = 0
		tot_fail = 0
		for index in region:
			tot += 1
			if index == -1:
				flag = False
				tot_fail += 1
				break
			else:
				x = vor.vertices[index, 0]
				y = vor.vertices[index, 1]
				z = vor.vertices[index, 2]
				# session.logger.info(str((x,y,z)))
				if not(bounding_box[0] - eps <= x and x <= bounding_box[1] + eps and
					   bounding_box[2] - eps <= y and y <= bounding_box[3] + eps and
					   bounding_box[4] - eps <= z and z <= bounding_box[5] + eps):
					flag = False
					tot_fail += 1
					# session.logger.info('failed')
					# break
		# session.logger.info('%d failed out of %d'%(tot_fail, tot))
		if flag:
			# session.logger.info('passed')
			regions.append(region)

			# find the point which lies inside
			points_to_filter.append(vor.points[vor.point_region == i][0,:])

	vor.filtered_points = np.array(points_to_filter)
	vor.filtered_regions = regions
	return vor

def centroid_region(vertices, map_2d, memory, threshold=0):
	global time_in_computation, time_point_computation, time_points_generation

	min_x = map_2d.shape[0] + 1
	min_y = map_2d.shape[1] + 1
	min_z = map_2d.shape[2] + 1
	max_x = -1
	max_y = -1
	max_z = -1
	for i in range(len(vertices)):
		min_x = min(min_x, vertices[i, 0])
		min_y = min(min_y, vertices[i, 1])
		min_z = min(min_z, vertices[i, 2])
		max_x = max(max_x, vertices[i, 0])
		max_y = max(max_y, vertices[i, 1])
		max_z = max(max_z, vertices[i, 2])
	min_x = int(min_x*map_2d.shape[0])
	max_x = int(max_x*map_2d.shape[0]) + 1
	min_y = int(min_y*map_2d.shape[1])
	max_y = int(max_y*map_2d.shape[1]) + 1
	min_z = int(min_z*map_2d.shape[2])
	max_z = int(max_z*map_2d.shape[2]) + 1
	max_x = min(max_x, map_2d.shape[0])
	max_y = min(max_y, map_2d.shape[1])
	max_z = min(max_z, map_2d.shape[2])
	min_x = max(min_x, 0)
	min_y = max(min_y, 0)
	min_z = max(min_z, 0)
	A = 0
	C_x = 0
	C_y = 0
	C_z = 0


	temp = np.zeros(vertices[0].shape)
	for i in range(len(vertices)):
		temp += vertices[i]/len(vertices)


	t = time.time()

	N_x = max_x - min_x
	N_y = max_y - min_y
	N_z = max_z - min_z
	grid = np.indices((N_x,N_y,N_z))
	g0 = (grid[0].ravel() + min_x + 0.5) /map_2d.shape[0]
	g1 = (grid[1].ravel() + min_y + 0.5) /map_2d.shape[1]
	g2 = (grid[2].ravel() + min_z + 0.5) /map_2d.shape[2]
	points = np.asarray([g0,g1,g2]).transpose()

	# points = []
	# index = 0
	# indices = np.zeros(map_2d.shape)
	# for i in range(min_x, max_x):
	# 	for j in range(min_y, max_y):
	# 		for k in range(min_z, max_z):
	# 			# points.append([(i + 0.5)/map_2d.shape[0], (j + 0.5)/map_2d.shape[1], (k + 0.5)/map_2d.shape[2]])
	# 			# if points1[index][0] != (i + 0.5)/map_2d.shape[0] or points1[index][1] != (j + 0.5)/map_2d.shape[1] or points1[index][2] != (k + 0.5)/map_2d.shape[2]:
	# 			# 	session.logger.info(str(N_x))
	# 			# 	session.logger.info(str(N_y))
	# 			# 	session.logger.info(str(N_z))
	# 			indices[i,j,k] = index
	# 			index += 1
	time_points_generation += time.time() -t
	# try:
	if len(points) > 0:
		t = time.time()
		in_array = Delaunay(vertices).find_simplex(points) >= 0
		time_in_computation += time.time() -t
	# except Exception as e:
		# raise Exception('ridam %s\n%s,\n%d,%d,%d,%d,%d,%d,%d,%d,%d'%(str(list(vertices)),str(list(points)),min_x,max_x,map_2d.shape[0],min_y,max_y,map_2d.shape[1],min_z,max_z,map_2d.shape[2],))
		t = time.time()

		map_array = map_2d[min_x:max_x, min_y:max_y, min_z:max_z].ravel()
		A += (in_array*map_array).sum()
		C_x += (in_array*map_array* g0).sum()
		C_y += (in_array*map_array* g1).sum()
		C_z += (in_array*map_array* g2).sum()

		# index = -1
		# for i in range(min_x, max_x):
		# 	for j in range(min_y, max_y):
		# 		for k in range(min_z, max_z):
		# 			index += 1
		# 			# if memory[i, j, k] == 1:
		# 			# 	continue
		# 			# if map_2d[i, j, k] < threshold:
		# 			# 	continue
		# 			# raise Exception('asdf %s'%(str(indices[i,j,k]),))
		# 			if in_array[index]:
		# 			# if in_array[int(indices[i,j,k])]:
		# 				# memory[i, j, k] = 1
		# 				A = A + map_2d[i,j,k]
		# 				C_x = C_x + map_2d[i,j,k] * (i + 0.5)/map_2d.shape[0]
		# 				C_y = C_y + map_2d[i,j,k] * (j + 0.5)/map_2d.shape[1]
		# 				C_z = C_z + map_2d[i,j,k] * (k + 0.5)/map_2d.shape[2]
		time_point_computation += time.time() - t
#     print('finished')
	if A == 0 :
		# if not temp.has_z:
		# 	raise Exception('ridi %s'%(str(list(vertices)),))
		return np.array([[temp[0], temp[1], temp[2]]]), A
	C_x /= A
	C_y /= A
	C_z /= A
#     print(np.array([[C_x, C_y]]))
	return np.array([[C_x, C_y, C_z]]), A

def plot(r, map_2d, bounding_box):
	t = time.time()
	vor = voronoi(r, bounding_box)
	session.logger.info('voronoi time %f'%(time.time()-t,))
	#TODO
	# return vor,0


	t = time.time()
	centroids = []
	weights = []
	memory = np.zeros(map_2d.shape)
	for region in vor.filtered_regions:
		vertices = vor.vertices[region + [region[0]], :]
		centroid, w = centroid_region(vertices, map_2d, memory)
		centroids.append(list(centroid[0, :]))
		weights.append(w)
	centroids = np.asarray(centroids)

	session.logger.info('centroid time %f'%(time.time()-t,))
	session.logger.info('time_in_computation time %f'%(time_in_computation,))
	session.logger.info('time_point_computation time %f'%(time_point_computation,))
	session.logger.info('time_points_generation time %f'%(time_points_generation,))
	return centroids, weights

def update(rob,centroids):
	interim_x = np.asarray(centroids[:,0] - rob[:,0])
	interim_y = np.asarray(centroids[:,1] - rob[:,1])
	interim_z = np.asarray(centroids[:,2] - rob[:,2])
	magn = [np.linalg.norm(centroids[i,:] - rob[i,:]) for i in range(rob.shape[0])]
	x = np.copy(interim_x)
	x = np.asarray([interim_x[i]/magn[i] for i in range(interim_x.shape[0])])
	y = np.copy(interim_y)
	y = np.asarray([interim_y[i]/magn[i] for i in range(interim_y.shape[0])])
	z = np.copy(interim_z)
	z = np.asarray([interim_z[i]/magn[i] for i in range(interim_z.shape[0])])
	nor = [np.linalg.norm([x[i],y[i],z[i]]) for i in range(x.shape[0])]
	temp = np.copy(rob)
	temp[:,0] = [rob[i,0] + 1*interim_x[i] for i in range(rob.shape[0])]
	temp[:,1] = [rob[i,1] + 1*interim_y[i] for i in range(rob.shape[0])]
	temp[:,2] = [rob[i,2] + 1*interim_z[i] for i in range(rob.shape[0])]
	return np.asarray(temp)

def get_init(map_3d, M, random_seed=None):
	cube_length = max(map_3d.shape)
	cube_length += cube_length%2
	map_3d = np.pad(map_3d, ((0, cube_length - map_3d.shape[0]), (0, cube_length - map_3d.shape[1]), (0, cube_length - map_3d.shape[2])), 'minimum')
	assert np.unique(map_3d.shape).size == 1, 'map must be cube, not non-cubic rectangular parallelepiped'
	N = map_3d.shape[0]
	assert N%2 == 0, 'N must be even'

	map_3d /= map_3d.sum() # 3d map to probability density
	map_3d_flat = map_3d.flatten()
	map_3d_idx = np.arange(map_3d_flat.shape[0])
	if random_seed is not None: np.random.seed(seed=random_seed)

	# this scales with M (the number of chosen items), not map_3d_idx (the possibilities to choose from)
	samples_idx = np.random.choice(map_3d_idx,size=M,replace=True,p=map_3d_flat) # chosen voxel indeces
	coords_1d = np.arange(-N//2,N//2)
	xyz = coords_n_by_d(coords_1d,d=3)
	rm0 = xyz[samples_idx] # pick out coordinates that where chosen. note that this assumes map_3d_idx matches with rows of xyz

	robs = []
	for i in range(len(rm0)):
		rob = [0,0,0]
		rob[0] = rm0[i][0]/N + 1/2 + np.random.normal(0, 0.1)/N
		rob[1] = rm0[i][1]/N + 1/2 + np.random.normal(0, 0.1)/N
		rob[2] = rm0[i][2]/N + 1/2 + np.random.normal(0, 0.1)/N
		robs.append(rob)
	robs = np.asarray(robs)

	return robs, map_3d


def iterate(map_3d, r0, max_iter=5):
	robots = r0
	bounding_box = np.array([0., 1., 0., 1., 0., 1.]) 

	for i in range(max_iter):
		centroids, weights = plot(robots,map_3d, bounding_box)
		robots = update(robots,centroids)

	#TODO
	# return centroids
	return robots