#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
from model.base import BaseXGB
from model.XGB_utils import\
    model_fit, grid_search, gen_xgb_input, model_inference, \
    construct_model_params


class XGBWrapper(BaseXGB):

    def __init__(self,
                 features,
                 target,
                 task,
                 nb_boost=100,
                 search_range=None,
                 hyper_params=None):
        super(XGBWrapper).__init__()
        self.features = features
        self.target = target
        self.task = task
        self.nb_boost = nb_boost
        if search_range is None:
            # use default
            search_range = {
                'max_depth': [3, 4, 5, 6, 7],
                'subsample': [0.5, 0.6, 0.7, 0.8, 0.9, 1]
            }
            if task == 'site' or task == 'binary':
                search_range['threshold'] = [
                    1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5
                ]
        self.search_range = search_range
        self.model_params = construct_model_params(hyper_params, task)

    def train(self, train_df):
        xg_train = gen_xgb_input(train_df, self.features, self.target)
        return model_fit(xg_train, self.model_params, self.nb_boost)

    def tune(self, train_df, val_df):
        best_model, best_hyperparam = grid_search(train_df, val_df,
                                                  self.features, self.target,
                                                  self.model_params,
                                                  self.nb_boost, self.task,
                                                  self.search_range)
        return best_model, best_hyperparam

    def predict(self, model, df):
        pred = model_inference(model, df, self.features)
        return pred
