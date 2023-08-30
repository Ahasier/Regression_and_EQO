function [crossValidatedCoefficients, squaredError, MSE, indexOfMSE, R2OutSample, R2InSample] = computeCrossValidation(abundanceData, functionalOutput, estimatedCoefficients, trainingSets, testSets, regressionMethod, beta0)
% COMPUTECROSSVALIDATION Performs cross-validation on the provided regression coefficients.
%
% INPUTS:
%   abundanceData: Matrix of abundance data.
%   functionalOuput: Values of functional outputs from different samples.
%   estimatedCoefficients: Matrix of estimated regression coefficients.
%   testSet: Indices or logical array indicating which samples to use as the test set.
%   regressionMethod: String specifying the regression method ('OLS' or others).
%   settings: Structure containing additional parameters and settings.
%
% OUTPUTS:
%   crossValidatedCoefficients: Coefficients obtained after cross-validation.
%   squaredError: Vector of squared errors from the cross-validation.
%   MSE: Mean squared error from the cross-validation.
%   indexOfMSE: Index of the mean squared error in the squaredError vector.

% Handle OLS regression method
if strcmp(regressionMethod, 'OLS')
    [crossValidateThresholds, squaredError, idxb, R2OutSample, R2InSample] = handleCrossValidationForOLS(abundanceData, functionalOutput, estimatedCoefficients, trainingSets, testSets, beta0);
else
    % Handle other regression methods
    [crossValidateThresholds, squaredError, idxb, R2OutSample, R2InSample] = handleCrossValidationForOthers(abundanceData, functionalOutput, estimatedCoefficients, trainingSets, testSets, beta0);
end
% Compute the minimum mean squared error and its index
[MSE, indexOfMSE] = min(squaredError(:));

% Obtain the coefficients corresponding to the minimum mean squared error
crossValidatedCoefficients = calculatecrossValidatedCoefficients(estimatedCoefficients, crossValidateThresholds, indexOfMSE);
end

%% Helper functions
% Perform cross-validation for Ordinary Least Squares regression
function [crossValidateThresholds, squaredError, idxb, R2OutSample, R2InSample] = handleCrossValidationForOLS(abundanceData, functionalOutput, estimatedCoefficients, trainingSets, testSets, beta0)
bs = 0:beta0/100:2*beta0;
crossValidateThresholds = 0:0.1:2*beta0;
len0 = length(bs);
len1 = size(estimatedCoefficients, 2);
len2 = length(crossValidateThresholds);
% Initialize arrays to store results
squaredError = zeros(1, len1);
idxb = zeros(1, len1);
R2InSample = zeros(1, len1);
R2OutSample = zeros(1, len1);

% Iterate through all cross-validation thresholds
for l = 1:len2
    % Iterate through all permutations
    for m = 1:len1
        allSquaredErrors = zeros(len0, len1);
        allR2InSample = zeros(len0, len1);
        allR2OutSample = zeros(len0, len1);
        % Iterate through all weight parameters (b)
        for n = 1:len0
            b = bs(n);
            threshold = crossValidateThresholds(l);
            BinaryBeta = estimatedCoefficients >= threshold;
            cappedBeta = BinaryBeta.*estimatedCoefficients;
            
            trainingSet = trainingSets(:, m);
            testSet = testSets(:, m);
            
            [allR2InSample(n,m), allR2OutSample(n,m), allSquaredErrors(n,m)] = computeSquaredErrorAndRSquared(abundanceData, functionalOutput, trainingSet, testSet, cappedBeta, b);
        end
    end
    % Get the squared error, r^2 correlation coefficients in and out of
    % sample, and their index of optimal weight at different threshold values.
    [squaredError(l), idxb(l)] = min(mean(allSquaredErrors,2));
    R2InSample(l) = mean(allR2InSample(idxb(l),:));
    R2OutSample(l) = mean(allR2OutSample(idxb(l),:));
end
end

% Handle other regression methods for cross-validation
function [crossValidateThresholds, squaredError, idxb, R2OutSample, R2InSample] = handleCrossValidationForOthers(abundanceData, functionalOutput, estimatedCoefficients, trainingSets, testSets, beta0)
bs = 0:beta0/10:2*beta0;
crossValidateThresholds = 0:0.1:2*beta0;
len0 = length(bs);
len1 = size(testSets, 2);
len2 = size(testSets, 3);
len3 = length(crossValidateThresholds);
squaredError = zeros(len1, len2);
idxb = zeros(len1, len2);
R2InSample = zeros(len1, len2);
R2OutSample = zeros(len1, len2);

% Iterate through all cross-validation thresholds
for l = 1:len3
    % Iterate through all regresssion (lambda) parameters
    for k = 1:len2
        % Iterate through all permutations
        for m = 1:len1
            allSquaredErrors = zeros(len0, len1);
            allR2InSample = zeros(len0, len1);
            allR2OutSample = zeros(len0, len1);
            % Iterate through all weight parameters (b)
            for n = 1:len0
                b = bs(n);
                threshold = crossValidateThresholds(l);
                BinaryBeta = estimatedCoefficients(:, k) >= threshold;
                cappedBeta = BinaryBeta.*estimatedCoefficients(:, k);
                
                trainingSet = trainingSets(:, m,k);
                testSet = testSets(:, m,k);
                
                [allR2InSample(n,m), allR2OutSample(n,m), allSquaredErrors(n,m)] = computeSquaredErrorAndRSquared(abundanceData, functionalOutput, trainingSet, testSet, cappedBeta, b);
            end
        end
        [squaredError(l, k), idxb(l, k)] = min(mean(allSquaredErrors, 2));
        R2InSample(l, k) = mean(allR2InSample(idxb(l, k), :));
        R2OutSample(l, k) = mean(allR2OutSample(idxb(l, k), :));
    end
end
end

% Obtain the cross-Validated coefficients corresponding to the minimum mean squared error
function crossValidatedCoefficients = calculatecrossValidatedCoefficients(estimatedCoefficients, crossValidateThresholds, indexOfMSE)
len = length(crossValidateThresholds);
thresholdedEstimatedCoefficients = zeros([len, size(estimatedCoefficients)]);
for l = 1:len
    binaryCoefficients = estimatedCoefficients > crossValidateThresholds(l);
    thresholdedEstimatedCoefficients(l,:) = binaryCoefficients.*estimatedCoefficients;
end
crossValidatedCoefficients = thresholdedEstimatedCoefficients(indexOfMSE,:);
crossValidatedCoefficients = reshape(crossValidatedCoefficients, [size(estimatedCoefficients, 1), 1]);
end

% Compute the R2 values and Squared Error for given data and coefficients
function [R2InSample, R2OutSample, SquaredError] = computeSquaredErrorAndRSquared(abundanceData, functionalOutput, trainingSet, testSet, cappedBeta, b)
% Compute predicted values for the training set
y_pred_insample = computePredictedValues(abundanceData(trainingSet, :), cappedBeta, b);
% Compute R2 for the training set
R2InSample = computeRSquared(functionalOutput(trainingSet), y_pred_insample);

% Compute predicted values for the test set
y_pred_outsample = computePredictedValues(abundanceData(testSet, :), cappedBeta, b);
% Compute squared error for the test set
SquaredError = computeSquaredErrors(functionalOutput(testSet), y_pred_outsample);
% Compute R2 for the test set
R2OutSample = computeRSquared(functionalOutput(testSet), y_pred_outsample);
end

% Compute the predicted values for given data, coefficients, and weight
function y_pred = computePredictedValues(abundanceData, coefficients, weight)
% Calculate the predicted values by matrix multiplication
y_pred = weight * abundanceData * coefficients;
end

% Compute the squared errors between the actual and predicted values
function SE = computeSquaredErrors(y_actual, y_pred)
% Calculate the sum of squared differences between actual and predicted values
SE = sum((y_actual - y_pred).^2);
end

% Compute the R2 (coefficient of determination) value
function R2 = computeRSquared(y_actual, y_pred)
% Compute sum of squared residuals
SSR = sum((y_actual - y_pred).^2);
% Compute total sum of squares
TSS = sum((y_actual - mean(y_actual)).^2);
% Calculate R2 as 1 minus the ratio of SSR to TSS
R2 = 1 - (SSR / TSS);
end