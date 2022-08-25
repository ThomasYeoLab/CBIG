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
        df (class DataFrame): Data
        mean (class Series): Mean
        std (class Series): Std
    """
    df = (df - mean) / std

    return df


def ICV_norm(df, ROIs):
    """
    Normalize ROIs by ICV

    Args:
        df (class DataFrame): Data
        ROIs (list): ROIs
    """
    df[ROIs] = df[ROIs].div(df['EstimatedTotalIntraCranialVol'].values, axis=0)
    return df
