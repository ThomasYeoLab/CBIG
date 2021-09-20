# Control analysis of using different primary gradients for parameterization
* Four different primary gradients are added in this test: inter-subject FC variability, first principal component of gene expression, first principal decomposition of structural covariance, second FC principal gradient
* `scripts` this folder contains the scripts to generate the results of different combinations of primary gradients
    * `CBIG_pMFM_basic_functions.py` this file contains all the necessary functions used in the pMFM analysis, such as using MFM to generate simulated fMRI signals; generating simulated FC and FCD.

    * `CBIG_pMFM_step1_training_Funcvar.py` this file is the wrapper to perform parameter estimation process. The parameterization are based only on group level inter-subject FC variability
    * `CBIG_pMFM_step2_validation_Funcvar.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step1_training_Funcvar.py`.
    * `CBIG_pMFM_step3_test_Funcvar.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step2_validation_Funcvar.py`.
    * `CBIG_pMFM_step4_training_FuncvarGrad.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level first FC gradient and inter-subject FC variability
    * `CBIG_pMFM_step5_validation_FuncvarGrad.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step4_training_FuncvarGrad.py`.
    * `CBIG_pMFM_step6_test_FuncvarGrad.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step5_validation_FuncvarGrad.py`. 
    * `CBIG_pMFM_step7_training_FuncvarT1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level T1w/T2w and inter-subject FC variability
    * `CBIG_pMFM_step8_validation_FuncvarT1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step7_training_FuncvarT1T2.py`.
    * `CBIG_pMFM_step9_test_FuncvarT1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step8_validation_FuncvarT1T2.py`. 

    * `CBIG_pMFM_step10_training_Gene.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on first principal component of gene expression
    * `CBIG_pMFM_step11_validation_Gene.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step10_training_Gene.py`.
    * `CBIG_pMFM_step12_test_Gene.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step11_validation_Gene.py`.
    * `CBIG_pMFM_step13_training_GeneGrad.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level first FC gradient and first principal component of gene expression
    * `CBIG_pMFM_step14_validation_GeneGrad.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step13_training_GeneGrad.py`.
    * `CBIG_pMFM_step15_test_GeneGrad.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step14_validation_GeneGrad.py`. 
    * `CBIG_pMFM_step16_training_GeneT1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level T1w/T2w and first principal component of gene expression
    * `CBIG_pMFM_step17_validation_GeneT1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step16_training_GeneT1T2.py`.
    * `CBIG_pMFM_step18_test_GeneT1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step17_validation_GeneT1T2.py`.

    * `CBIG_pMFM_step19_training_Struct.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on first principal decomposition of structural covariance
    * `CBIG_pMFM_step20_validation_Struct.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step19_training_Struct.py`.
    * `CBIG_pMFM_step21_test_Struct.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step20_validation_Struct.py`.
    * `CBIG_pMFM_step22_training_StructGrad.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level first FC gradient and first principal decomposition of structural covariance
    * `CBIG_pMFM_step23_validation_StructGrad.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step22_training_StructGrad.py`.
    * `CBIG_pMFM_step24_test_StructGrad.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step23_validation_StructGrad.py`. 
    * `CBIG_pMFM_step25_training_StructT1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level T1w/T2w and first principal decomposition of structural covariance
    * `CBIG_pMFM_step26_validation_StructT1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step25_training_StructT1T2.py`.
    * `CBIG_pMFM_step27_test_StructT1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step26_validation_StructT1T2.py`. 

    * `CBIG_pMFM_step28_training_GradPC2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on second FC principal gradient
    * `CBIG_pMFM_step29_validation_GradPC2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step28_training_GradPC2.py`.
    * `CBIG_pMFM_step30_test_GradPC2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step29_validation_GradPC2.py`.
    * `CBIG_pMFM_step31_training_GradPC2Grad.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level first FC gradient and second FC principal gradient
    * `CBIG_pMFM_step32_validation_GradPC2Grad.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step31_training_GradPC2Grad.py`.
    * `CBIG_pMFM_step33_test_GradPC2Grad.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step32_validation_GradPC2Grad.py`. 
    * `CBIG_pMFM_step34_training_GradPC2T1T2.py` this file is the wrapper to perform parameter estimation process. The parameterization are based on group level T1w/T2w and second FC principal gradient
    * `CBIG_pMFM_step35_validation_GradPC2T1T2.py` this file is the wrapper to perform estimated model parameter validation process based on results of `CBIG_pMFM_step34_training_GradPC2T1T2.py`.
    * `CBIG_pMFM_step36_test_GradPC2T1T2.py` this file is the wrapper to perform estimated model parameter test process based on results of `CBIG_pMFM_step35_validation_GradPC2T1T2.py`. 

# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part2_pMFM_control_analysis/Primary_gradients/scripts` and run the scripts.
* For python scripts, please first activate the environment (e.g. `source activate pMFM`) and then directly run the scripts.
