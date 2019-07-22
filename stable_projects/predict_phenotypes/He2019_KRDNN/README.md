## Reference

He, T., Kong, R., Holmes, A., Nguyen, M., Sabuncu, M., Eickhoff, S.B., Bzdok, D., Feng, J. and Yeo, B.T., 2019. [**Deep Neural Networks and Kernel Regression Achieve Comparable Accuracies for Functional Connectivity Prediction of Behavior and Demographics**](https://www.biorxiv.org/content/10.1101/473603v1), under review.

----
## Background

There is significant interest in the development and application of deep neural networks (DNNs) to neuroimaging data. Here, we compared the performance of three DNN architectures and a classical machine learning algorithm (kernel regression) in predicting individual phenotypes from whole-brain resting-state functional connectivity (RSFC) patterns.

----

## Code Release
### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_He2019_KRDNN](https://github.com/ThomasYeoLab/Standalone_He2019_KRDNN)

### Download whole repository
If you want to use the code from our lab's other stable projects (other than He2019_KRDNN), you would need to download the whole CBIG repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.13.2-He2019_KRDNN](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.13.2-He2019_KRDNN)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b He2019_KRDNN v0.13.2-He2019_KRDNN
    ```

### Usage
- Our code uses MATLAB and Python, here are info about MATLAB and Python setup:
	- MATLAB: we tested our code in MATLAB r2018b and r2014a (the example results for kernel regression cross validation are a little bit different between r2014a and r2018b. Other than this, all the example and replication results are same for r2014a and r2018b)
	- Python
		1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/distribution/#download-section) with Python 3.x if you don't have conda
		2. Create conda environment from our `config/CBIG_KRDNN_python_env.yml` file by `conda env create -f config/CBIG_KRDNN_python_env.yml`
		3. If you have keras installed, replace yours `$HOME/.keras/keras.json` with `config/keras.json`. If not, create `$HOME/.keras/keras.json` with content of `config/keras.json`
- The example of our code is detailed in `examples/README.md`
- If you have access to HCP and UK Biobank dataset, you can replicate our result using the instructions detailed in `replication/README.md`.


----

## Updates
- Release v0.13.0 (10/07/2019): Initial release of He2019_KRDNN project
- Release v0.13.2 (22/07/2019):
    1. Update RNG generator in MATLAB for compatibility
    2. Update title of paper in README files
    3. Add information about hyperparameter in replication/README.md

----

## Bugs and Questions

Please contact He Tong at hetong1115@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

