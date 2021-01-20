#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', '-f', nargs='+', required=True)
    parser.add_argument('--suffix', '-s', required=True)

    return parser.parse_args()


def combine_frames(frames):
    for idx in range(1, len(frames)):
        assert np.all(np.asarray(frames[0]['RID']) == np.asarray(frames[idx]['RID']))
        assert np.all(np.asarray(frames[0]['Forecast Month']) == np.asarray(frames[idx]['Forecast Month']))
        assert np.all(np.asarray(frames[0]['Forecast Date']) == np.asarray(frames[idx]['Forecast Date']))

    out_frame = frames[0].copy()
    out_frame['Forecast Date'] = out_frame['Forecast Date'].apply(lambda x: '-'.join(x.split('-')[:2]))
    out_frame['CN relative probability'] = np.mean(np.array([f['CN relative probability'] for f in frames]), axis=0)
    out_frame['MCI relative probability'] = np.mean(np.array([f['MCI relative probability'] for f in frames]), axis=0)
    out_frame['AD relative probability'] = np.mean(np.array([f['AD relative probability'] for f in frames]), axis=0)

    adas = np.mean(np.array([f.ADAS13 for f in frames]), axis=0)
    out_frame.ADAS13 = adas
    out_frame['ADAS13 50% CI lower'] = np.clip(adas - 6, 0, 85)
    out_frame['ADAS13 50% CI upper'] = np.clip(adas + .5, 0, 85)

    norm_vent = np.mean(np.array([f.Ventricles_ICV for f in frames]), axis=0)
    out_frame.Ventricles_ICV = norm_vent
    out_frame['Ventricles_ICV 50% CI lower'] = np.clip(norm_vent - 2e-3, 0., 1.)
    out_frame['Ventricles_ICV 50% CI upper'] = np.clip(norm_vent + 2e-3, 0., 1.)

    return out_frame


if __name__ == '__main__':
    args = get_args()

    final_frame = combine_frames([pd.read_csv(ff) for ff in args.files])
    final_frame.to_csv('TADPOLE_Submission_%s.csv' % args.suffix, index=False)
