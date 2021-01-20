#!/usr/bin/env python
from __future__ import print_function
import csv
import os
import threading


class OptimBase(object):
    def __init__(self, executable, fixed_args, log_file):
        self.executable = executable
        self.fixed_args = fixed_args
        self.log_file = log_file
        assert not os.path.isfile(self.log_file), '%s exits. Choose a different name' % self.log_file
        self.best_result = 100
        self.lock = threading.Lock()
        self.f_eval_count = 0

        self.hyper_map = None
        self.dim = 0
        self.func = []

        self.init_param_list()

        with open(self.log_file, 'w') as fhandler:
            header = ['bestScore', 'score', 'nb_eval', 'fname'] + self.hyper_map.keys()
            csv.writer(fhandler).writerow(header)

    def init_param_list(self):
        raise NotImplementedError()

    def param_vec2dict(self, params):
        return {k: self.func[v](params[v]) for k, v in self.hyper_map.items()}

    def log(self, result, nb_eval, fname, var_args):
        row = [self.best_result, result, nb_eval, fname]
        row += [var_args[k] for k in self.hyper_map.keys()]
        with open(self.log_file, 'a') as fhandler:
            csv.writer(fhandler).writerow(row)

    def execute(self, var_args):
        raise NotImplementedError()
