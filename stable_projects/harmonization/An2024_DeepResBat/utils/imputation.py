#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np


def forward_filling(df, feature):
    """
    Forward filling
    """
    mean = np.nanmean(df[feature].values)
    if feature == 'DX':
        mean = 0.
    subs = np.unique(df.RID)
    for sub in subs:
        sub_mask = df.RID == sub
        dates = df[sub_mask]['EXAMDATE'].values
        dates = np.sort(dates)
        bl_date = dates[0]
        bl_mask = (sub_mask) & (df.EXAMDATE == bl_date)
        bl_value = df.loc[bl_mask, [feature]].values[0][0]
        if np.isnan(bl_value):
            bl_value = mean
            df.loc[bl_mask, [feature]] = mean
        ff_value = bl_value  # value to forward fill
        for date in dates[1:]:
            date_mask = (sub_mask) & (df.EXAMDATE == date)
            date_value = df.loc[date_mask, [feature]].values[0][0]
            if np.isnan(date_value):
                df.loc[date_mask, [feature]] = ff_value
            else:
                # update ff_value
                ff_value = date_value
    return df


def ff_multi_cols(df, cols):
    """
    Run forward filling for multiple columns on d dataframe
    """
    for col in cols:
        df = forward_filling(df, col)
    return df
