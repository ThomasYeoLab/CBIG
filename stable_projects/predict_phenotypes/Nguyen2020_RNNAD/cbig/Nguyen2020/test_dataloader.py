# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
import unittest
import numpy as np
import cbig.Nguyen2020.dataloader as dataloader


class DataloaderTest(unittest.TestCase):
    """ Unit tests for different filling strategies """

    def setUp(self):
        self.xin = np.array(
            [0, 6, 36, 60, 66, 72, 78, 84, 90, 96, 102, 109, 120],
            dtype=np.int)
        self.yin = np.array([
            np.nan, np.nan, np.nan, 8., np.nan, 7., np.nan, 8., np.nan, 5.,
            np.nan, np.nan, 7.
        ])
        self.xout = np.hstack([np.arange(0, 60, 30), np.arange(60, 144, 6)])

    def test_bl_fill(self):
        y_bl_fill, _, _ = dataloader.bl_fill(self.xin, self.yin, 0., self.xout)
        ref = np.zeros(len(y_bl_fill), dtype=bool)
        ref[[1, 3, 5, 7, 9, 10, 11, 13, 14, 15]] = True
        assert np.all(np.isnan(y_bl_fill) == ref)

    def test_ff_fill(self):
        y_ff_fill, _, _ = dataloader.ff_fill(self.xin, self.yin, 0., self.xout)
        assert np.any(~np.isnan(y_ff_fill))

    def test_ln_fill_partial(self):
        y_ln_fill_partial, _, _ = dataloader.ln_fill_partial(
            self.xin, self.yin, 0., self.xout)
        ref = np.zeros(len(y_ln_fill_partial), dtype=bool)
        ref[[1, 10, 11, 13, 14, 15]] = True
        assert np.all(np.isnan(y_ln_fill_partial) == ref)

    def test_ln_fill_full(self):
        y_ln_fill_full, _, _ = dataloader.ln_fill_full(self.xin, self.yin, 0.,
                                                       self.xout)
        ref = np.zeros(len(y_ln_fill_full), dtype=bool)
        ref[[1, 13, 14, 15]] = True
        assert np.all(np.isnan(y_ln_fill_full) == ref)
