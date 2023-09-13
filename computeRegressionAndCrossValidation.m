function [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% COMPUTEREGRESSIONANDCROSSVALIDATION integrates the regression and cross-validation processes.
% 
% INPUTS:
%   abundanceData: Matrix containing the abundance data of different taxa across samples.
%   functionalOutput: Vector containing the specific functional output in the samples.
%   numPermutations: Number of times to perform the integrated process.
%   regressionMethod: String specifying the regression method ('L0', 'OLS', etc.).
%   settings: A structure with additional parameters and settings.
%
% OUTPUTS:
%   Various metrics related to the integrated regression and cross-validation process, including coefficients, out-of-sample errors, R-squared values, and optimal thresholds.

% Check if cross-validating over different threshold
if isfield(settings, 'Threshold') && strcmp(settings.Threshold, 'cv')
    [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationViaThresholding(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);
    
elseif isnan(settings.Threshold) % Check for AIC Setting
    [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationViaAIC(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);
    
else % Otherwise, using a given threshold
    [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationGivenThreshold(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);
end
end