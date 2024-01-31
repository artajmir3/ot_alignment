import math
import numpy as np
import random
import time
import json

import mrcfile
import pandas as pd

import torch

from scipy.spatial.transform import Rotation

from Bio import PDB
from Bio.PDB.vectors import Vector, rotmat

from .point_sampler import *
from .transform import *
from .ot import *

from chimerax.core.commands import run






def build_coup(x, y, z, x1, y1, z1, pi, verbose=False):
	res = np.array(pi)
	all_coup = []
	for i in range(len(x1)):
		maxi = None
		maxx = 0
		s = 0
		for j in range(len(x)):
			s += res[j,i]
			if res[j,i] > maxx:
	#             print(maxx, res[i,j])
	#             print(i,j)
				maxx = res[j,i]
				maxi = j
	#     print('!')
		if verbose:
			print(i,maxi, maxx, s)
	#     if maxi/maxx >0.8:
		if maxi is not None:
	#     if i == 0 and maxi >1e-4:
			all_coup.append((i, maxi))
	return all_coup

def find_optimal_alignment(x, y, z, x1, y1, z1, all_coup, verbose=False):
	A = []
	B = []
	for i in range(len(all_coup)):
		X = [x[all_coup[i][1]], y[all_coup[i][1]], z[all_coup[i][1]]]
		X = np.array(X)
		A.append(X)
	#     print(all_coup[i][1])
		X = [x1[all_coup[i][0]], y1[all_coup[i][0]], z1[all_coup[i][0]]]
		X = np.array(X)
		B.append(X)
	Abar = np.array([0., 0., 0.])
	Bbar = np.array([0., 0., 0.])
	for i in range(len(all_coup)):
		Abar += 1/len(all_coup) * A[i]
		Bbar += 1/len(all_coup) * B[i]
	if verbose:
		print(Abar)
		print(Bbar)
	H = np.zeros((3,3))
	for i in range(len(all_coup)):
		H += np.outer(A[i] - Abar, (B[i] - Bbar))
	# print(H)
	U, S, V = np.linalg.svd(H)
	# print(U,S,V)
	# R = np.matmul(np.transpose(V), np.transpose(U))
	R = np.matmul(U,V)
	d = np.linalg.det(R)
#     R = np.matmul(np.matmul(U,np.diag([1,1,d])),V)
	if verbose:
		print(H)
		print(np.diag(S))
		print(np.matmul(np.matmul(U,np.diag(S)),V))
		print(V)
		print(R)
		print(np.linalg.det(R))
	r = Rotation.from_matrix(R)
	return Abar, Bbar, r

def plot_alignment(x, y, z, x1, y1, z1, Abar, Bbar, r, all_coup):
	x2 = []
	y2 = []
	z2 = []
	x3 = []
	y3 = []
	z3 = []
	for i in range(len(x)):
		x2.append(x[i] - Abar[0])
		y2.append(y[i] - Abar[1])
		z2.append(z[i] - Abar[2])
		x3.append(x1[i] - Bbar[0])
		y3.append(y1[i] - Bbar[1])
		z3.append(z1[i] - Bbar[2])
		
	x3, y3, z3 = perform(x3, y3, z3, change_quat_format(r.as_quat()))
	
	plt.rcParams["figure.figsize"] = (20,20)
	fig = plt.figure()
	# ax = fig.gca(projection='3d', adjustable='box')
	ax = fig.add_subplot(projection='3d')
	ax.scatter(x3, y3, z3,  marker='o')
	ax.scatter(x2, y2, z2,  marker='o')
	for i in range(len(all_coup)):
		ax.plot([x2[all_coup[i][1]], x3[all_coup[i][0]]], [y2[all_coup[i][1]], y3[all_coup[i][0]]], zs=[z2[all_coup[i][1]], z3[all_coup[i][0]]], c='k', alpha=0.05)
	ax.set_xlim(-60,60)
	ax.set_ylim(-60,60)
	ax.set_zlim(-60,60)
	plt.show()



def align_pointcloud(pointcloud1, pointcloud2, num=2, verbose=False):
	x, y, z = pointcloud1.get_global_coords()
	x1, y1, z1 = pointcloud2.get_global_coords()
	alignments = []
	for i in range(num):
		if verbose:
			print(i)
			
#         if verbose:
#             fig = plt.figure()
#             # ax = fig.gca(projection='3d', adjustable='box')
#             ax = fig.add_subplot(projection='3d')
#             ax.scatter(x1, y1, z1,  marker='o')
#             ax.scatter(x, y, z,  marker='o', alpha=0.5)    
#             plt.show()
			
		t_perm  = None
		if num > 1:
			perm = np.random.permutation(len(x))
			t_perm = np.zeros((len(x), len(x)))
			for i in range(len(x)):
				t_perm[i, perm[i]] = 1/len(x)
			t_perm = torch.from_numpy(t_perm)
		
		pi, gamma, dist = ugw(x, y, z, x1, y1, z1, init=t_perm)
		
		all_coup = build_coup(x, y, z, x1, y1, z1, pi)
		
		Abar, Bbar, r = find_optimal_alignment(x, y, z, x1, y1, z1, all_coup)
		
		if verbose:
			plot_alignment(x, y, z, x1, y1, z1, Abar, Bbar, r, all_coup)
		
		covering = [0] * len(x)
		for i in range(len(all_coup)):
			covering[all_coup[i][1]] += 1
		
		alignments.append({'Abar':Abar, 'Bbar':Bbar, 'r':r, 'score':dist, 'covering':covering})
	return alignments

def find_best_score(alignments):
	best_score = float('inf')
	best = None
	for alignment in alignments:
		if alignment['score'] < best_score:
			best = alignment
			best_score = alignment['score']
	return best
