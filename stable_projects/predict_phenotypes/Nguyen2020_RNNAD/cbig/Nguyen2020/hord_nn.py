#!/usr/bin/env python
from __future__ import print_function
import subprocess32 as subprocess

import numpy as np

from hord_base import OptimBase


def dict2paramlist(mapping):
    return list(sum([('--%s' % k, v) for k, v in mapping.items()], tuple()))


class NetOptim(OptimBase):
    def execute(self, var_args):
        args = self.executable + self.fixed_args
        args += dict2paramlist(var_args)
        args = [p for p in args if len(p) > 0]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        out, err = proc.communicate()
        out = out.strip()
        assert proc.returncode == 0, err
        assert len(out.split('\n')) == 1, ('Expect only one line', out)
        pred_fname, mauc, bca, adas, vent = out.split(',')
        score = self._calc_score(float(mauc), float(bca), float(adas), float(vent))

        return pred_fname, score

    def objfunction(self, params, gpu_id):
        assert len(params) == self.dim, 'Dimension mismatch'

        var_args = self.param_vec2dict(params)
        var_args['gpu'] = str(gpu_id)

        pred_fname, result = self.execute(var_args)

        with self.lock:
            if self.best_result > result:
                self.best_result = result
                print('Params:', *dict2paramlist(var_args))
                print('Score:', result)

            self.f_eval_count += 1
            nb_eval = self.f_eval_count

        self.log(result, nb_eval, pred_fname, var_args)

        return result


def param_dropout():
    return 0, 5, lambda x: '%.1f' % (x / 10.)


def param_learning_rate():
    return -5, -2, lambda x: '%e' % (10**x)


def param_l2_decay():
    return -7, -4, lambda x: '%e' % (10**x)


def param_hid_dim():
    return 7, 9, lambda x: '%d' % (2**x)


def param_no_layers():
    return 1, 3, lambda x: '%d' % x


class tadpole_objective(NetOptim):
    def init_param_list(self):
        self.dim = 6
        self.continuous = np.arange(0, 1)
        self.integer = np.arange(1, self.dim)

        # Hyperparameter to optimise:
        params = ['lr', 'weight_decay', 'h_drop', 'i_drop', 'h_size', 'nb_layers']
        self.hyper_map = {p: i for i, p in enumerate(params)}
        m = self.hyper_map

        self.xlow = np.zeros(self.dim)
        self.xup = np.zeros(self.dim)
        self.func = [None for _ in range(self.dim)]

        index = m['lr']
        self.xlow[index], self.xup[index], self.func[index] = param_learning_rate()

        index = m['weight_decay']
        self.xlow[index], self.xup[index], self.func[index] = param_l2_decay()

        index = m['h_drop']
        self.xlow[index], self.xup[index], self.func[index] = param_dropout()

        index = m['i_drop']
        self.xlow[index], self.xup[index], self.func[index] = param_dropout()

        index = m['h_size']
        self.xlow[index], self.xup[index], self.func[index] = param_hid_dim()

        index = m['nb_layers']
        self.xlow[index], self.xup[index], self.func[index] = param_no_layers()

    def _calc_score(self, mauc, bca, adas, vent):  # lower is better
        total = -mauc - bca + adas + vent
        print('\t', total, mauc, bca, adas, vent)
        return total if np.isfinite(total) else 100000
