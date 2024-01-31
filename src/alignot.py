import math
import numpy as np
import random
import time
import json


from .point_sampler import *
from .transform import *
from .ot import *

	
def dot(a, b):
	s = 0
	for i in range(4):
		s += a[i] * b[i]
	return s

	
def alignot(point_cloud1, point_clouds2, lr=0.005, max_iter=100, reg=0.1, num_samples=1, verbose=False, random_seed=None):
	# return [[1,0,0,0]], [0]
	# max_iter = 2
	x, y, z = point_cloud1.get_global_centered_coords() 
	xr, yr, zr = point_clouds2.get_global_centered_coords()
	px = Quaternion(Polynomial([Term(coef=0, exps=[0,0,0,0,1,0,0,0])]), Polynomial([Term(coef=1, exps=[0,0,0,0,0,1,0,0])]), 
				Polynomial([Term(coef=1, exps=[0,0,0,0,0,0,1,0])]), Polynomial([Term(coef=1, exps=[0,0,0,0,0,0,0,1])]))

	bx = q*px*qs
	i_der = []
	j_der = []
	k_der = []
	for i in range(4):
		i_der.append(bx.i_pol.derivative(i))
		j_der.append(bx.j_pol.derivative(i))
		k_der.append(bx.k_pol.derivative(i))

	if random_seed is not None:
		random.seed(random_seed)
	
	vals = get_quaternion_vals(0, 0, 0, 1)
	quaternions = []
	costs = []
	OT_time = 0
	grad_time = 0
	sample_time = 0
	rotate_time = 0
	prev = None
	u = None
	v = None
	for i in range(max_iter):
		t = time.time()
		if verbose:
			if i % 10 == 9:
				print('Iteration number %d, the wasserstein deistance is %.2f'%(i, costs[-1]))
		quaternions.append([])
		for j in range(4):
			quaternions[-1].append(vals[j])
		
		xx,yy,zz = perform(x, y, z, vals)
		rotate_time += time.time() - t
		
		t = time.time()
		if u is not None:
			prev = (u,v)
		T,cost,u,v = OT(xr,yr,zr,xx,yy,zz,reg=reg,prev=prev)
		costs.append(cost)
		OT_time += time.time() -t
		t = time.time()
		
		for s in range(num_samples):
			sample_point = int(random.random()* len(x))
			
			

			x1 = x[sample_point]
			y1 = y[sample_point]
			z1 = z[sample_point]
			dest = None
			maxx = 0
			for j in range(len(x)):
				if T[sample_point][j] * len(x) > maxx:
					maxx = T[sample_point][j] * len(x)
					dest = j
			x2 = xr[dest]
			y2 = yr[dest]
			z2 = zr[dest]
			
			
			norm_grad = 0
			sample_time += time.time() - t
			t = time.time()

			grad = [0, 0, 0, 0]

			new_vals = []
			for k in range(4):
				new_vals.append(vals[k])
			new_vals.append(0)
			new_vals.append(x1)
			new_vals.append(y1)
			new_vals.append(z1)

			i_vals = bx.i_pol.evaluate(new_vals)
			i_der_vals = []
			j_vals = bx.j_pol.evaluate(new_vals)
			j_der_vals = []
			k_vals = bx.k_pol.evaluate(new_vals)
			k_der_vals = []
			for k in range(4):
				i_der_vals.append(i_der[k].evaluate(new_vals))
				j_der_vals.append(j_der[k].evaluate(new_vals))
				k_der_vals.append(k_der[k].evaluate(new_vals))
			

			for j in range(len(x)):
				if j != dest:
					continue
				x2 = xr[j]
				y2 = yr[j]
				z2 = zr[j]

				temp = [0, 0, 0, 0]
				for k in range(4):

					temp[k] += 2 * i_der_vals[k] * (i_vals - x2)
					temp[k] += 2 * j_der_vals[k] * (j_vals - y2)
					temp[k] += 2 * k_der_vals[k] * (k_vals - z2)

					grad[k] += temp[k] * T[sample_point][j] * len(x) /num_samples
					
		d_prod = dot(grad, vals)
		for j in range(4):
			grad[j] -= d_prod * vals[j]
		
		for j in range(4):
			vals[j] -= lr * grad[j]
		
		norm = math.sqrt(vals[0]**2 + vals[1]**2 + vals[2]**2 + vals[3]**2)
		for j in range(4):
			vals[j] /= norm
		
		norm_grad = math.sqrt(grad[0]**2 + grad[1]**2 + grad[2]**2 + grad[3]**2)
		grad_time += time.time() -t

	if verbose:
		print('Time spent for optimal transport is ' + str(OT_time) + ' second(s).')
		print('Time spent for computing gradient is ' + str(grad_time) + ' second(s).')
		print('Time spent for rotating is ' + str(rotate_time) + ' second(s).')
		print('Time spent for sampling is ' + str(sample_time) + ' second(s).')
		print('Final cost: ' + str(costs[-1]))
	return quaternions, costs, [OT_time, grad_time, rotate_time, sample_time]




