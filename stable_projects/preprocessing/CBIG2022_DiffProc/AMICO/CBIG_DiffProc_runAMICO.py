#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Leon Ooi and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import sys
import os
import amico


def runAMICO(output_dir, subj, dwi_dir, mask):
    # setup
    amico.core.setup()
    ae = amico.Evaluation("AMICO", subj)

    # create scheme file
    bvals = os.path.join(dwi_dir, subj, subj + ".bval")
    bvecs = os.path.join(dwi_dir, subj, subj + ".bvec")
    if not os.path.isdir(os.path.join(output_dir, 'scheme_files')):
        os.mkdir(os.path.join(output_dir, 'scheme_files'))
    diffusion_scheme = os.path.join(output_dir, 'scheme_files',
                                    subj + ".scheme")
    amico.util.fsl2scheme(bvals, bvecs, schemeFilename=diffusion_scheme)

    # load data
    dwi_file = os.path.join(dwi_dir, subj, subj + ".nii.gz")
    if not os.path.isfile(dwi_file):
        print("Could not fine nii.gz file, trying nii format")
        dwi_file = os.path.join(dwi_dir, subj, subj + ".nii")
        print(dwi_file)
        if not os.path.isfile(dwi_file):
            print(
                "ERROR: DWI file is not in nii.gz nor nii format, exiting...")
            exit()
    ae.load_data(
        dwi_filename=dwi_file,
        scheme_filename=diffusion_scheme,
        mask_filename=mask,
        b0_thr=0)

    # generate kernels
    os.chdir(os.path.join(output_dir, 'output'))
    ae.set_model("NODDI")
    ae.generate_kernels()
    ae.load_kernels()

    # generate results
    ae.fit()
    ae.save_results()


if __name__ == '__main__':
    runAMICO(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
