"""
Test packages' installation.
"""

import pytest

class TestClass:
    def test_mapca(self):
        import mapca

    def test_nibabel(self):
        import nibabel

    def test_nilearn(self):
        import nilearn

    def test_numpy(self):
        import numpy

    def test_sklearn(self):
        import sklearn    

    def test_scipy(self):
        import scipy
    
    def test_tedana(self):
        import tedana

    def test_torch(self):
        import torch

print("Everything seems fine. Your Python environment is set up properly!")
