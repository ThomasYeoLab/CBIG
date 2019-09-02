# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
import unittest
from datetime import datetime
from dateutil.relativedelta import relativedelta

import numpy as np

import cbig.Nguyen2020.misc as misc


def inefficient_month_between(end, start):
    """ Find duration (in months) between two dates """
    series = np.array([relativedelta(months=i) for i in range(240)])
    series += start

    return np.abs(series - end).argmin()


class MiscTest(unittest.TestCase):
    """ Unit test for finding duration (in months) between two dates """

    def test_month_between(self):
        start = datetime.strptime('2010-01-12', '%Y-%m-%d')

        for i in range(31):
            end = datetime.strptime('2010-03-%2.d' % (i + 1), '%Y-%m-%d')
            assert misc.month_between(end, start) == inefficient_month_between(
                end, start)
