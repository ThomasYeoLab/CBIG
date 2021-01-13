"""
Test packages' installation.
"""

import pytest

class TestClass:
    def test_nibabel(self):
        import nibabel as nib

    def test_nilearn(self):
        import nilearn

    def test_pydicom(self):
        import dicom

    def test_nipy(self):
        import nipy

print("Everything seems fine. Your Python environment is set up properly!")
