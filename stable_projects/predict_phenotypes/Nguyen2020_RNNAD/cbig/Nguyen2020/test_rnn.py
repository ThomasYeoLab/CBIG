# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
import unittest
import torch
import cbig.Nguyen2020.rnn as rnn


class RnnCellTest(unittest.TestCase):
    """ Unit tests for recurrent cells """

    def setUp(self):
        torch.manual_seed(0)
        self.in_features = 10
        self.hidden_size = 20
        self.batch_size = 3
        self.length = 15

    def test_MinimalRNNCell(self):
        cell = rnn.MinimalRNNCell(self.in_features, self.hidden_size)
        seq = torch.randn(self.length, self.batch_size, self.in_features)
        h_t = torch.randn(self.batch_size, self.hidden_size)
        for i in range(self.length):
            h_t = cell(seq[i], h_t)

        self.assertAlmostEqual(h_t.sum().item(), -3.026607, 6)

    def test_LssCell(self):
        cell = rnn.LssCell(self.in_features, self.hidden_size)
        seq = torch.randn(self.length, self.batch_size, self.in_features)
        h_t = torch.randn(self.batch_size, self.hidden_size)
        for i in range(self.length):
            h_t = cell(seq[i], h_t)

        self.assertAlmostEqual(h_t.sum().item(), 60.245380, 6)
