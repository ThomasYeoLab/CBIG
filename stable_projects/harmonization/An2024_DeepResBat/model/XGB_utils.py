#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import copy
import itertools
import xgboost as xgb
from sklearn.metrics import accuracy_score, mean_squared_error
from utils.metrics import data_pred_metric
import warnings

warnings.filterwarnings('ignore')


def gen_xgb_input(df, features, target=None):
    """
    Generate input file for XGBoost models
    """
    if target is None:
        x = df[features].values
        return xgb.DMatrix(x, feature_names=features)

    else:
        x, y = df[features].values, df[target].values
        return xgb.DMatrix(x, y, feature_names=features)


def model_fit(xg_train, model_params, nb_boosts=100):
    """
    Fit a XGBoost model & make prediction on validation set
    """
    if 'threshold' in model_params.keys():
        xgb_model_params = copy.deepcopy(model_params)
        del xgb_model_params['threshold']
        # fit model
        model = xgb.train(params=xgb_model_params,
                          dtrain=xg_train,
                          num_boost_round=nb_boosts,
                          verbose_eval=False)
    else:
        # fit model
        model = xgb.train(params=model_params,
                          dtrain=xg_train,
                          num_boost_round=nb_boosts,
                          verbose_eval=False)

    return model


def model_inference(model, df, features):
    """
    Inference using trained model
    """
    return model.predict(gen_xgb_input(df, features))


def prob2value(pred_prob, threshold=None):
    """
    Transfer probility prediction to real values
    """
    if threshold is None:
        # multi-class
        assert pred_prob.shape[1] > 2, "Binary pred prob needs threshold!"
        return pred_prob.argmax(axis=1)
    else:
        # binary
        return (pred_prob >= threshold) * 1


def grid_search(train_df, val_df, features, target, params, nb_boosts, task,
                search_range):
    """
    Perform grid search for XGBoost models
    """
    # generate input for XGBoost model
    xg_train = gen_xgb_input(train_df, features, target)
    y_val = val_df[target].values

    # get initial performance with default hyper-parameters
    best_model = model_fit(xg_train, params, nb_boosts)
    val_pred = model_inference(best_model, val_df, features)
    if task == 'binary':
        best_score = metric_func(prob2value(val_pred, params['threshold']),
                                 y_val, task)
    elif task == 'site':
        _, best_score = data_pred_metric(val_df['RID'].values, val_pred, y_val,
                                         params['threshold'])
    elif task == 'multi':
        best_score = metric_func(prob2value(val_pred), y_val, task)
    else:
        best_score = metric_func(val_pred, y_val, task)

    # grid search for hyper-parameters
    def get_hyper_params_grid(search_range):
        """
        Get a grid combination for hyper_pramas
        """
        keys, values = zip(*search_range.items())
        grids = [dict(zip(keys, v)) for v in itertools.product(*values)]

        return grids, keys

    # note that threshold is a hyper-params no need for train
    if task == 'site' or task == 'binary':
        search_range_ = copy.deepcopy(search_range)
        del search_range['threshold']
    hyper_parmas_combs, hyper_parmas = get_hyper_params_grid(search_range)
    # initialize best hyper-parameters
    best_hyperparam = dict()
    for hyper_parma in hyper_parmas:
        best_hyperparam[hyper_parma] = params[hyper_parma]
    best_hyperparam['threshold'] = 0.01 # A starting point for best threshold

    # grid search
    for hyper_params_comb in hyper_parmas_combs:
        for hyper_parma in hyper_parmas:
            params[hyper_parma] = hyper_params_comb[hyper_parma]
        # get prediction
        model = model_fit(xg_train, params, nb_boosts)
        val_pred = model_inference(model, val_df, features)
        if task == 'binary' or task == 'site':
            best_threshold_score = 0
            best_threshold = search_range_['threshold'][0]
            for threshold in search_range_['threshold']:
                hyper_params_comb['threshold'] = threshold
                if task == 'binary':
                    score = metric_func(
                        prob2value(val_pred, hyper_params_comb['threshold']),
                        y_val, task)
                else:
                    _, score = data_pred_metric(val_df['RID'].values, val_pred,
                                                y_val,
                                                hyper_params_comb['threshold'])
                if score > best_threshold_score:
                    best_threshold_score = score
                    best_threshold = threshold
            if best_threshold_score > best_score:
                best_score = best_threshold_score
                best_model = model
                best_hyperparam = copy.deepcopy(hyper_params_comb)
                best_hyperparam['threshold'] = best_threshold
        else:
            if task == 'multi':
                score = metric_func(prob2value(val_pred), y_val, task)
            else:
                score = metric_func(val_pred, y_val, task)
            if score > best_score:
                best_score = score
                best_model = model
                best_hyperparam = copy.deepcopy(hyper_params_comb)

    return best_model, best_hyperparam


def metric_func(pred, gt, task='regress'):
    """
    Get metric for given prediction and ground_truth, larger <=> better
    """
    if task == 'regress':
        return -1 * mean_squared_error(pred.reshape(
            (-1, 1)), gt.reshape((-1, 1)))
    else:
        return accuracy_score(pred.reshape((-1, 1)), gt.reshape((-1, 1)))


def construct_model_params(hyper_params, task='regress'):
    """
    Construct model params by feeding hyper_params
    """
    if hyper_params is None:
        model_params = {
            'booster': 'gbtree',
            'nthread': 1,
            'max_depth': 6,
            'subsample': 1.0,
        }
        if task == 'regress':
            model_params['objective'] = 'reg:squarederror'
            model_params['eval_metric'] = 'rmse'
        elif task == 'multi':
            model_params['objective'] = 'multi:softprob'
            model_params['eval_metric'] = 'merror'
            model_params['num_class'] = 3
        else:
            model_params['objective'] = 'binary:logistic'
            model_params['eval_metric'] = 'error'
            model_params['threshold'] = 0.01

        return model_params
    else:
        # with hyper-parameters
        model_params = {'booster': 'gbtree', 'nthread': 1}
        for key in list(hyper_params.keys()):
            model_params[key] = hyper_params[key]
        if task == 'regress':
            model_params['objective'] = 'reg:squarederror'
            model_params['eval_metric'] = 'rmse'
        elif task == 'multi':
            model_params['objective'] = 'multi:softprob'
            model_params['eval_metric'] = 'merror'
            model_params['num_class'] = 3
        else:
            model_params['objective'] = 'binary:logistic'
            model_params['eval_metric'] = 'error'

        return model_params
