## Reference

+ Nguyen, M., Sun, N., Alexander D.C., Feng J., Yeo B.T.T., 2018. [**Modeling Alzheimer’s disease progression using deep recurrent neural networks**](https://doi.org/10.1109/prni.2018.8423955), PRNI, 2018.
+ Nguyen, M., He T., An L., Alexander D.C., Feng J., Yeo B.T.T., 2020. [**Predicting Alzheimer’s disease progression using deep recurrent neural networks**](https://doi.org/10.1016/j.neuroimage.2020.117203), NeuroImage, 117203.

----
## Background

Early identification of people at risk of developing Alzheimer’s disease (AD) would be beneficial for developing treatments.
This project uses recurrent neural network (RNN) to predict the progression of AD in subjects over the long term.
Temporal interpolation strategies are used to deal with missing data, thus making efficient use of longitudinal data.

----

## Code Release
### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_Nguyen2020_RNNAD](https://github.com/ThomasYeoLab/Standalone_Nguyen2020_RNNAD)

### Download whole repository
If you want to use the code from our lab's other stable projects (other than Nguyen2020_RNNAD), you would need to download the whole CBIG repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.22.2-Update_He2022_MM](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.22.2-Update_He2022_MM)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b Nguyen2020_RNNAD v0.22.2-Update_He2022_MM
    ```

### Usage
- This code is compatible with Python 2.7, to create a Python environment similar to what was used for this project:
    1. Install [Anaconda](https://www.anaconda.com/distribution/#download-section) with Python 2.7
    2. Create Anaconda environment from our `replication/config/CBIG_RNN_python_env.yml` file by `conda env create -f replication/config/CBIG_RNN_python_env.yml -n NAME_OF_ENVIRONMENT`. This should create an environment compatible to the latest version of the code.
    3. However, if you want to replicate the experiments in the papers, use the `replication/config/CBIG_RNN_python_env.yml` configuration file instead.
- An example of how to use the code is detailed in `examples/README.md`.


----

## Updates
- Release v0.22.2 (16/03/2022): 
    1. Fix the typo in examples/CBIG_RNN_example.sh
    2. Fix the typo in tadpole/CBIG_ensemble_prediction.sh
- Release v0.18.1 (20/01/2021):
    1. Update unit test to accommodate to the new HPC.
    2. Release scripts to run hyperparameter search for Nguyen2020_RNNAD project.
- Release v0.14.0 (02/09/2019): Initial release of Nguyen2020_RNNAD project

----

## Bugs and Questions

Please contact Minh Nguyen at minhnb3192@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
