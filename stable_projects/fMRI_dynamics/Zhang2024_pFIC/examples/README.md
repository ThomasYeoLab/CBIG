## Examples

This folder contains demo input data, configuration files and reference outputs. The reference outputs are generated using **NVIDIA RTX3090 with CUDA version 11.7**. To exactly replicate the results, please make sure that you are using the same hardwares and same versions of Python packages. The details of the Python environment can be found under `model/pFIC.yml`. 


### Usage
To run the demo script, run
```
cd script
CBIG_pFIC_wrapper_example.sh
```
The script first activates the dedicated pFIC Python environment, then it calls the pFIC model training scripts and saves out the training results under the `examples/output` folder, which will only be created upon running this script.

This script serves as both an example and a unit-test test case, thus it generates the outputs and automatically compares them with the reference outputs. 

Typical run-time using a single RTX3090 card is about 30 minutes.