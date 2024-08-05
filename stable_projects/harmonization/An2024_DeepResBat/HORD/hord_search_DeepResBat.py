#!/usr/bin/env python
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

from __future__ import print_function
import argparse
import sys
import numpy as np
from poap.controller import ThreadController
from poap.controller import BasicWorkerThread
from HORD.pySOTpy2.experimental_design import LatinHypercube
from HORD.pySOTpy2.rbf_interpolant import RBFInterpolant
from HORD.pySOTpy2.rbf_surfaces import CubicRBFSurface
from HORD.pySOTpy2.sot_sync_strategies import SyncStrategyNoConstraints

from HORD.pySOT_ext import BoundaryDYCORS
from HORD.hord_nn_DeepResBat import tadpole_objective


def main(args):
    np.random.seed(0)
    data = tadpole_objective(args.executable.split(), sys.argv[1:],
                             args.log_file)

    nb_candidates = args.multiplier * data.dim
    smpl_method = BoundaryDYCORS(data=data, numcand=nb_candidates)
    exp_design = LatinHypercube(dim=data.dim, npts=(data.dim + 1))

    print('Number of threads:', args.nthreads)
    print('Maximum number of evaluations:', args.maxeval)
    print('No. candidates:', nb_candidates)
    print('Search strategy:', smpl_method.__class__.__name__)
    print('Experimental design:', exp_design.__class__.__name__)
    print('Surrogate: Cubic RBF')

    # Create a strategy and a controller
    controller = ThreadController()
    surface = RBFInterpolant(surftype=CubicRBFSurface, maxp=args.maxeval)
    controller.strategy = \
        SyncStrategyNoConstraints(
            worker_id=0,
            data=data,
            maxeval=args.maxeval,
            nsamples=args.nthreads,
            exp_design=exp_design,
            response_surface=surface,
            sampling_method=smpl_method
        )

    # Launch the threads and give them access to the objective function
    def create_function(n):
        def function(x):
            return data.objfunction(x, gpu_id=n)
        return function
    for i in range(args.nthreads):
        n = str(i % args.nGPUs)
        function = create_function(n)
        controller.launch_worker(BasicWorkerThread(controller, function))

    # Run the optimization strategy
    result = controller.run()

    print('Best value found:', result.value)
    print('Best solution found:', data.param_vec2dict(result.params[0]))


def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=str, default='')
    parser.add_argument('--nthreads', type=int, default=1)
    parser.add_argument('--nGPUs', type=int, default=4)
    parser.add_argument('--maxeval', type=int, default=200)
    parser.add_argument('--log_file', type=str, default='')
    parser.add_argument('--multiplier', type=int, default=100)
    return parser


if __name__ == '__main__':
    args, _ = get_arg_parser().parse_known_args()
    main(args)
