"""
Test packages' installation.
"""

import pytest

class TestClass:
    def test_nibabel(self):
        import nibabel as nib

    def test_nilearn(self):
        import nilearn

    def test_nipype(self):
        import nipype.interfaces.fsl as fsl
        import nipype.pipeline.engine as pe
        from nipype import Node, JoinNode, Workflow

    def test_pydicom(self):
        import dicom

    def test_nipy(self):
        import nipy

print("Everything seems fine. Your Python environment is set up properly!")
