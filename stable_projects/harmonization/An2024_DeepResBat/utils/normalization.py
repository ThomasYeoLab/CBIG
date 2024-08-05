#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''


def z_norm(df, mean, std):
    """
    z normalization

    Args:
        df (_type_): _description_
        mean (_type_): _description_
        std (_type_): _description_

    Returns:
        _type_: _description_
    """
    df = (df - mean) / std

    return df


def ICV_norm(df, ROIs):
    """
    Normalize ROIs by ICV

    Args:
        df (_type_): _description_
        ROIs (_type_): _description_

    Returns:
        _type_: _description_
    """
    df[ROIs] = df[ROIs].div(df['EstimatedTotalIntraCranialVol'].values, axis=0)
    return df


def min_max_norm(df, min_value, max_vale):
    """
    min-max normalization

    Args:
        df (_type_): _description_
        min_value (_type_): _description_
        max_vale (_type_): _description_

    Returns:
        _type_: _description_
    """
    df = (df - min_value) / (max_vale - min_value)

    return df
