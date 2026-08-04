"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch


def check_predictor_results():
    reference_output_path = './reference_output/predictor/predictor_1/predictor_record/predictor_record.pth'
    output_path = './output/predictor/predictor_1/predictor_record/predictor_record.pth'
    reference_output = torch.load(reference_output_path, weights_only=False)
    output = torch.load(output_path, weights_only=False)

    difference = torch.sum(torch.abs(output['mse_mean'] - reference_output['mse_mean']))
    if difference <= 1e-2:
        print('Predictor results are the same as the reference outputs.')
        return 0
    else:
        print('Predictor results differ from the reference outputs.')
        return 1


if __name__ == "__main__":
    predictor_result = check_predictor_results()
    if predictor_result == 0:
        print('Passed! All results are the same as the reference outputs.')
    else:
        raise Exception('Failed. Results differ too much from the reference outputs.')
