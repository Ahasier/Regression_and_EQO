
# README for Regression_and_EQO

## Project Description
This MATLAB project uses regression methods to recover microbial functional groups (coefficient == 1 if some taxa is included in that functional group, and coefficient == 0 if not) from mock or real datasets. The data includes an abundance table of different taxa in different samples, and the specific functional output in different samples. it recovers the coefficients relating to whether-or-not a taxa belongs to that functional group and evaluates the prediction accuracy of various regression models under different conditions.

## Directory Structure and Main Files

```
Regression_and_EQO/
│
├── computeCrossValidation.m                      # Script for cross-validation
├── computeRegression.m                           # Script for regression models
├── computeRegressionAndCrossValidation.m         # Unified script for regression and cross-validation
├── generateData.m                                # Script for data generation
├── saveResults.m                                 # Script for saving results
├── computeAndSaveRegressionResults.m             # Script for computing and saving regression results
├── calculateAccuracy.m                           # Script for calculating accuracy metrics
├── loadRealData.m                                # Script for loading real dataset
├── solveUnifiedRegression.m                      # Script for unified regression model
├── generateSyntheticData.m                       # Script for generating synthetic data
├── formulateOptimization.m                       # Script for formulating optimization problems
├── RunMe.m                                       # Main script to run the entire analysis
├── solveOLSRegression.m                          # Script for Ordinary Least Squares (OLS) regression
│
├── Meta_Tara.csv                                 # Real dataset in CSV format
├── tree100taxa.mat                               # MATLAB data file for synthetic phylogenetic tree
├── tree100taxaReal.mat                           # MATLAB data file for real phylogenetic tree
│
├── results/                                      # Directory for storing result files
│   ├── Betas/                                    # Subdirectory for regression coefficient files
│   └── temporaryStorage/                         # Subdirectory for temporary storage
│       └── Betas/                                # Subdirectory for temporary regression coefficient files
│
└── functions/                                    # Directory for utility and helper functions
    ├── computeRSquared.m
    ├── computeSquaredError.m
    ├── getCompleteAccuracyFilePath.m
    ├── handleExtraPhylogeneticFeatures.m
    ├── loadOrInitializeAccuracyResultsFile.m
    ├── plotGapMap.m
    ├── plotMap.m
    ├── plotstyle.m
    ├── setOptionsAndNames.m
    ├── usingRealData.m
    └── visualizeResults.m
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

## Configuration settings
`configurations/config.json` stores the configuration settings for Regression & EQO. The default settings are:
```
{
    "Beta0": 1, // "Beta0" is the average non-zero value when generating ground truth coefficients.
    "BetaEps": 0, // "BetaEps" is the variance of non-zero value when generating ground truth.
    "Threshold": "NaN", // "Threshold" is determines the binarization method. If Threshold = nan, use AIC method; if Threshold is a given value, use that value to threshold recovered coefficients; else if threshold is 'cv', use cross-validation to determine the best threshold value. If using EQO, no need for this field.
    "requirePositivity": "Off", // Whether require positivity in regressions or not. If using EQO, no need for this field.
    "RealAbd": "Off", // Whether generating mimic-real abundance data or not.
    "weight": "On", // Whether weight estimated coefficients in computeAIC or not. if not doing regressions with AIC method, no need for this field.
    "DiagnosticMod": "Off", // If "DiagnosticMod" is "On", storing extra results for diagnostics.
    "maxLambda": 10, // Maximum number of lambda parameter in the cross-validation of LASSO/ Ridge/ other regressions. If using OLS or EQO, no need for this field.
    "usePhylogeny": "Off" // If "usePhylogeny" is "On", incorporating phylogenetic information to regressions. Not needed if using EQO, or if generated synthetic data is phylogenetically-independent (phylogenyDependency = 0).
}
```