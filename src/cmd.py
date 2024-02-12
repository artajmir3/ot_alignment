# vim: set expandtab shiftwidth=4 softtabstop=4:
import time

from .point_sampler import *
from .empot import *
from .alignot import *


from chimerax.core.commands import CmdDesc
from chimerax.map import MapsArg, MapStepArg, Float1or3Arg, ValueTypeArg
from chimerax.core.commands import IntArg, StringArg, FloatArg, BoolArg
from chimerax.core.errors import UserError as CommandError

from chimerax.core.commands import run
from chimerax.geometry import Place



def perform_empot(session, volumes, n=500, thresh=0, num=2, random_seed=None, sampling_method='trn', local_refinement=True):
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
	v_res = transform_map(volumes[1], best_alignment)

	if local_refinement:
		run(session, 'fitmap #%d inMap #%d'%(v_res.id[0], volumes[0].id[0]))

	session.logger.info("spent %.2f(s)"%(time.time() - t,))

def perform_alignot(session, volumes, n=500, thresh=0, lr=0.0001, max_iter=500, reg=100, random_seed=None, sampling_method='trn', local_refinement=True):
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
	v_res = transform_map(volumes[1], alignment)

	if local_refinement:
		run(session, 'fitmap #%d inMap #%d'%(v_res.id[0], volumes[0].id[0]))

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



varg = [('volumes', MapsArg)]
empot_desc = CmdDesc(required = varg, 
					keyword = [('n', IntArg),
							('thresh', FloatArg),
							('num', IntArg),
							('random_seed', IntArg),
							('sampling_method', StringArg),
							('local_refinement', BoolArg)], 
					synopsis = 'Perform EMPOT to align given volumes')
alignot_desc = CmdDesc(required = varg, 
					keyword = [('n', IntArg),
							('thresh', FloatArg),
							('max_iter', IntArg),
							('random_seed', IntArg),
							('lr', FloatArg),
							('reg', FloatArg),
							('sampling_method', StringArg),
							('local_refinement', BoolArg)], 
					synopsis = 'Perform AlignOT to align given volumes')
illustrate_points_desc = CmdDesc(required = varg, 
							keyword = [('n', IntArg),
									('size', IntArg),
									('random_seed', IntArg),
									('thresh', FloatArg),
									('color', StringArg),
									('sampling_method', StringArg)], 
							synopsis = 'Generate and illustrate point clouds')