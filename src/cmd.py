# vim: set expandtab shiftwidth=4 softtabstop=4:
import time

from .point_sampler import *
from .empot import *
from .alignot import *


from chimerax.core.commands import CmdDesc
from chimerax.map import MapsArg, MapStepArg, Float1or3Arg, ValueTypeArg
from chimerax.core.commands import IntArg, StringArg, FloatArg
from chimerax.core.errors import UserError as CommandError

from chimerax.core.commands import run
from chimerax.geometry import Place



def hello_world(session, volumes):
# 	# All command functions are invoked with ``session`` as its
# 	# first argument.  Useful session attributes include:
# 	#   logger: chimerax.core.logger.Logger instance
# 	#   models: chimerax.core.models.Models instance
	session.logger.info("Hello world!")
# 	session.logger.info(str(volumes))
	session.logger.info(str(volumes[0]))
	session.logger.info(str(volumes[0].scene_position))
	from chimerax.geometry import identity
	point_to_scene_transform = identity()
	session.logger.info(str(help(identity())))
	m, xyz_to_ijk_tf = volumes[0].matrix_and_transform(point_to_scene_transform, subregion = None, step = None)
	session.logger.info(str(m))
	session.logger.info(str(xyz_to_ijk_tf))
	session.logger.info(str(help(xyz_to_ijk_tf)))

	trn = TRNPointSampler(volumes[0], random_seed=1)
	trn.map_th = m.copy()
	trn.sample(50)
	points = []
	for i in range(len(trn.x)):
		points.append(np.asarray([trn.x[i], trn.y[i], trn.z[i]]))
	points = np.asarray(points)
	xyz_to_ijk_tf.inverse().transform_points(points, in_place = True)
	trn.x = []
	trn.y = []
	trn.z = []
	for i in range(len(points)):
		trn.x.append(points[i][0])
		trn.y.append(points[i][1])
		trn.z.append(points[i][2])
	size=2
	color='#ff5a5a'
	for i in range(len(points)):
		run(session, "shape sphere radius %d center %f,%f,%f color %s"%(size, trn.x[i], trn.y[i], trn.z[i], color))


	v = volumes[0].writable_copy(require_copy = True, copy_colors = False)
	p = Place(axes=np.asarray([[1,0,0],[0,1,0],[0,0,1]]), origin=np.asarray([10,0,0]))
	v.scene_position = p * v.scene_position


# 	session.logger.info(str(volumes[1]))
# 	session.logger.info(str(help(volumes[0])))
# 	session.logger.info(str(help(volumes[0].session)))
# 	#session.logger.info(str(volumes[0].data.rotation))
# 	#session.logger.info(str(volumes[0].data.origin))
# 	#session.logger.info(str(volumes[0].data.voxel_volume()))
# 	#session.logger.info(str(volumes[0].data.ijk_to_xyz([0,0,0])))
# 	#session.logger.info(str(volumes[0].ijk_to_global_xyz([0,0,0])))
# 	#session.logger.info(str(volumes[0].ijk_to_global_xyz([0.5,0,0])))
# 	#session.logger.info(str(volumes[0].data.ijk_to_xyz([50,50,50])))
# 	#session.logger.info(str(volumes[0].data.matrix()[50][50][50]))
# 	#session.logger.info(str(volumes[0].matrix()[50][50][50]))
# 	#session.logger.info(str(volumes[0].id))
# 	session.logger.info("Hello world!")

# 	trn = TRNPointSampler(volumes[0])
# 	trn.threshold(0.557)
# 	trn.sample(500)
# 	#trn.show_points()

# 	trn1 = TRNPointSampler(volumes[1])
# 	trn1.threshold(0.557)
# 	trn1.sample(500)
# 	#trn1.show_points()

# 	alignments = align_pointcloud(trn, trn1)

# 	best_alignment = find_best_score(alignments)
# 	session.logger.info(str(best_alignment))
# 	session.logger.info("Hello world!")
# 	transform_map(volumes[1], best_alignment)

def perform_empot(session, volumes, n=500, thresh=0, num=2, random_seed=None, sampling_method='trn'):
	session.logger.info("perform empot")
	t = time.time()


	if len(volumes) != 2:
		raise CommandError('this command requires exactly 1 volume, got %d' % len(volumes))

	sampling_class = None
	if sampling_method == 'trn':
		sampling_class = TRNPointSampler
	elif sampling_method == 'cvt':
		sampling_class = CVTPointSampler
	else:
		raise CommandError('sampling method should be either trn or cvt, got %d' % len(volumes))

	alignments = []
	for i in range(num):
		if random_seed is not None:
			random_seed += 1

		trn = sampling_class(volumes[0], random_seed=random_seed)
		trn.threshold(thresh)
		trn.sample(n)
		trn1 = sampling_class(volumes[1], random_seed=random_seed)
		trn1.threshold(thresh)
		trn1.sample(n)

		alignments += align_pointcloud(trn, trn1, num=2)	

	best_alignment = find_best_score(alignments)
	session.logger.info(str(best_alignment))
	transform_map(volumes[1], best_alignment)

	session.logger.info("spent %.2f(s)"%(time.time() - t,))

def perform_alignot(session, volumes, n=500, thresh=0, lr=0.0001, max_iter=500, reg=100, random_seed=None, sampling_method='trn'):
# def perform_alignot(session, volumes, n=500, thresh=0, lr=0.000005, max_iter=100, reg=100000, random_seed=None):
	session.logger.info("perform alignot")
	t = time.time()

	if len(volumes) != 2:
		raise CommandError('this command requires exactly 1 volume, got %d' % len(volumes))

	sampling_class = None
	if sampling_method == 'trn':
		sampling_class = TRNPointSampler
	elif sampling_method == 'cvt':
		sampling_class = CVTPointSampler
	else:
		raise CommandError('sampling method should be either trn or cvt, got %d' % len(volumes))
		
	trn = sampling_class(volumes[0], random_seed=random_seed)
	trn.threshold(thresh)
	trn.sample(n)
	trn1 = TRNPointSampler(volumes[1], random_seed=random_seed)
	trn1.threshold(thresh)
	trn1.sample(n)

	quaternions, costs, times = alignot(trn, trn1, lr=lr, max_iter=max_iter, reg=reg, random_seed=random_seed)

	session.logger.info(str(costs))
	session.logger.info(str(times))

	alignment = {}
	alignment['Abar'] = trn.get_centroid()
	alignment['Bbar'] = trn1.get_centroid()
	alignment['q'] = quaternions[-1]
	alignment['q'][0] *= -1
	session.logger.info(str(alignment))
	transform_map(volumes[1], alignment)
	# transform_map1(volumes[1], alignment)

	session.logger.info("spent %.2f(s)"%(time.time() - t,))

def illustrate_points(session, volumes, n=500, random_seed=None, size=2, color='#ff5a5a', thresh=0, sampling_method='trn'):
	session.logger.info("illustrate points")
	t = time.time()

	if len(volumes) != 1:
		raise CommandError('this command requires exactly 1 volume, got %d' % len(volumes))

	trn = None
	if sampling_method == 'trn':
		trn = TRNPointSampler(volumes[0], random_seed=random_seed)
	elif sampling_method == 'cvt':
		trn = CVTPointSampler(volumes[0], random_seed=random_seed)
	else:
		raise CommandError('sampling method should be either trn or cvt, got %d' % len(volumes))
	trn.threshold(thresh)
	trn.sample(n)
	trn.show_points(size=size, color=color)

	session.logger.info("spent %.2f(s)"%(time.time() - t,))



# CmdDesc contains the command description.  For the
# "hello" command, we expect no arguments.
varg = [('volumes', MapsArg)]
hello_world_desc = CmdDesc(required = varg)
empot_desc = CmdDesc(required = varg, 
					keyword = [('n', IntArg),
							('thresh', FloatArg),
							('num', IntArg),
							('random_seed', IntArg),
							('sampling_method', StringArg)], 
					synopsis = 'Perform EMPOT to align given volumes')
alignot_desc = CmdDesc(required = varg, 
					keyword = [('n', IntArg),
							('thresh', FloatArg),
							('max_iter', IntArg),
							('random_seed', IntArg),
							('lr', FloatArg),
							('reg', FloatArg),
							('sampling_method', StringArg)], 
					synopsis = 'Perform AlignOT to align given volumes')
illustrate_points_desc = CmdDesc(required = varg, 
							keyword = [('n', IntArg),
									('size', IntArg),
									('random_seed', IntArg),
									('thresh', FloatArg),
									('color', StringArg),
									('sampling_method', StringArg)], 
							synopsis = 'Generate and illustrate point clouds')