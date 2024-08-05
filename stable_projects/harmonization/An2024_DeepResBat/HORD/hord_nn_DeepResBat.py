#!/usr/bin/env python
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

from __future__ import print_function
import subprocess32 as subprocess
import numpy as np

from HORD.hord_base import OptimBase
from config import global_config


def dict2paramlist(mapping):
    return list(sum([('--%s' % k, v) for k, v in mapping.items()], tuple()))


class NetOptim(OptimBase):

    def execute(self, var_args):
        args = self.executable + self.fixed_args
        args += dict2paramlist(var_args)
        args = [p for p in args if len(p) > 0]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        out, err = proc.communicate()
        out = out.decode('utf-8')
        out = out.strip()
        assert proc.returncode == 0, err
        assert len(out.split('\n')) == 1, ('Expect only one line', out)
        valROIMSE, valSiteAcc = out.split(',')
        score = self._calc_score(float(valROIMSE), float(valSiteAcc))

        return score

    def objfunction(self, params, gpu_id):
        assert len(params) == self.dim, 'Dimension mismatch'

        var_args = self.param_vec2dict(params)
        var_args['GPU'] = str(gpu_id)

        result = self.execute(var_args)

        with self.lock:
            if self.best_result > result:
                self.best_result = result
                print('Params:', *dict2paramlist(var_args))
                print('Score:', result)

            self.f_eval_count += 1
            nb_eval = self.f_eval_count

        self.log(result, nb_eval, var_args)

        return result


def param_learning_rate():
    return -2, -1, lambda x: '%e' % (10**x)


def param_dropout():
    return 0., 0.5, lambda x: '%.1f' % (x)


def param_alpha():
    return 0.01, 1., lambda x: '%.2f' % (x)


def param_lam():
    return 0.01, 1., lambda x: '%.3f' % (x)


def param_gam():
    return 0.01, 10., lambda x: '%.2f' % (x)


def param_lrstep():
    return 10, 999, lambda x: '%d' % (x)


def param_latentdim():
    return 5, 9, lambda x: '%d' % (2**x)


def param_h1():
    return 5, 9, lambda x: '%d' % (2**x)


def param_h2():
    return 5, 9, lambda x: '%d' % (2**x)


def param_h3():
    return 5, 9, lambda x: '%d' % (2**x)


def param_h4():
    return 5, 9, lambda x: '%d' % (2**x)


def param_nb_layers():
    return 2, 4, lambda x: '%d' % x


class tadpole_objective(NetOptim):
    """
    Class to find best hyperparameters point on given search space
    """

    def init_param_list(self):
        self.dim = global_config.hord_dim
        self.continuous = np.arange(0, global_config.hord_continous_dim)
        self.integer = np.arange(global_config.hord_continous_dim, self.dim)

        # Hyperparameter to optimise:
        params = [
            'lr', 'drop_out', 'alpha', 'lambda_', 'gamma', 'lr_step',
            'latent_dim', 'h1', 'h2', 'h3', 'h4', 'nb_layers'
        ]
        self.hyper_map = {p: i for i, p in enumerate(params)}
        m = self.hyper_map

        self.xlow = np.zeros(self.dim)
        self.xup = np.zeros(self.dim)
        self.func = [None for _ in range(self.dim)]

        index = m['lr']
        self.xlow[index], self.xup[index], self.func[index] = \
            param_learning_rate()

        index = m['drop_out']
        self.xlow[index], self.xup[index], self.func[index] = \
            param_dropout()

        index = m['alpha']
        self.xlow[index], self.xup[index], self.func[index] = \
            param_alpha()

        index = m['lambda_']
        self.xlow[index], self.xup[index], self.func[index] = \
            param_lam()

        index = m['gamma']
        self.xlow[index], self.xup[index], self.func[index] = \
            param_gam()

        index = m['lr_step']
        self.xlow[index], self.xup[index], self.func[index] = param_lrstep()

        index = m['latent_dim']
        self.xlow[index], self.xup[index], self.func[index] = param_latentdim()

        index = m['h1']
        self.xlow[index], self.xup[index], self.func[index] = param_h1()

        index = m['h2']
        self.xlow[index], self.xup[index], self.func[index] = param_h2()

        index = m['h3']
        self.xlow[index], self.xup[index], self.func[index] = param_h3()

        index = m['h4']
        self.xlow[index], self.xup[index], self.func[index] = param_h4()

        index = m['nb_layers']
        self.xlow[index], self.xup[index], self.func[index] = param_nb_layers()

    def _calc_score(self, valROIMSE, valSiteAcc):  # lower is better
        total = valROIMSE / 2 + valSiteAcc
        print('\t', total, valROIMSE, valSiteAcc)
        return total if np.isfinite(total) else 100000
