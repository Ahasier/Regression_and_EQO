
# README for Regression_and_EQO

## Project Description
This MATLAB project uses regression methods to recover functional groups (coefficient == 1 if some taxa is included in that functional group, and coefficient == 0 if not) from mock or real datasets. The data includes an abundance table of different taxa in different samples, and the specific functional output in different samples. it recovers the coefficients relating to whether-or-not a taxa belongs to that functional group and evaluates the prediction accuracy of various regression models under different conditions.

## Directory Structure and Main Files

```
Regression_and_EQO/
│
├── computeCrossValidation.m            # Script for cross-validation
├── generateData.m                      # Script for data generation
├── saveResults.m                       # Script for saving results
├── computeAndSaveRegressionResults.m   # Script for computing and saving regression results
├── calculateAccuracy.m                 # Script for calculating accuracy metrics
├── loadRealData.m                      # Script for loading real dataset
├── solveUnifiedRegression.m            # Script for unified regression model
├── generateSyntheticData.m             # Script for generating synthetic data
├── formulateOptimization.m             # Script for formulating optimization problems
├── RunMe.m                             # Main script to run the entire analysis
├── solveOLSRegression.m                # Script for Ordinary Least Squares (OLS) regression
├── computeRegression.m                 # Script for computing regression models
│
├── Meta_Tara.csv                       # Real dataset in CSV format
├── tree100taxa.mat                     # MATLAB data file related to synthetic phylogenetic tree
├── tree100taxaReal.mat                 # MATLAB data file  related to real phylogenetic tree
│
├── results/                            # Directory for storing result files
│   └── Betas/                          # Subdirectory for storing regression coefficient files
└── functions/                          # Directory for utility and helper functions
    ├── handleExtraPhylogeneticFeatures.m
    ├── generatePermutations.m
    ├── indicesOfAccuracyMatrixElement.m
    ├── loadOrInitializeAccuracyResultsFile.m
    ├── getCompleteAccuracyFilePath.m
    ├── usingRealData.m
    ├── setOptionsAndNames.m
    ├── plotMap.m
    ├── plotGapMap.m
    ├── plotstyle.m
```

## How to Run the Project
1. Open MATLAB and navigate to the `Regression_and_EQO` directory.
2. Run the `RunMe.m` script. This script is presumed to be the main entry point for the project:
   ```
   run('RunMe.m')
   ```
3. The `RunMe.m` script is expected to call the necessary functions and scripts to perform the analysis, and save the results in the appropriate directories.

## Dependencies
- MATLAB (The version used for this project is not specified, so it is recommended to use a recent version of MATLAB.)
