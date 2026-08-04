"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch

if __name__ == "__main__":
    reference_output_path = './reference_output/pfic/test/seed1/test_results.pth'
    output_path = './output/pfic/test/seed1/test_results.pth'
    reference_output = torch.load(reference_output_path, weights_only=False)
    output = torch.load(output_path, weights_only=False)

    parameter_difference = torch.sum(torch.abs(output['parameter'] - reference_output['parameter']))
    output_test_loss = output['corr_loss'] + output['l1_loss'] + output['ks_loss']
    reference_output_test_loss = reference_output['corr_loss'] + reference_output['l1_loss'] + reference_output[
        'ks_loss']
    loss_difference = torch.sum(torch.abs(output_test_loss - reference_output_test_loss))
    if parameter_difference + loss_difference <= 1e-6:
        print('Passed! Results are the same as the reference outputs.')
    else:
        raise Exception('Failed. Results differ too much from the reference outputs.')
