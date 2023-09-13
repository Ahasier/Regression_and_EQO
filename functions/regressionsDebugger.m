function [results, coefficientsComparison] = regressionsDebugger(regressionMethod, numTaxa, numSamples)
% Load result file
path = SetPathsForDataAndResults('betaResults');
numPermutations = 100;
varargin = {'Beta0', 1, 'BetaEps', 0.5, 'Threshold', 'cv', 'RealAbd','On','DiagnosticMod', 'On', 'requirePositivity', 'On'};
[settings, fullIdentifier] = setOptionsAndNames(varargin{:});
resultfile = [path, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numTaxa), '_nSpl', num2str(numSamples), '.mat'];
load(resultfile);

% If cross-validate via thresholding,check the out of sample R^2 and MSE
if isfield(settings, 'Threshold') && strcmp(settings.Threshold, 'cv')
    for n = 1:numPermutations
        optimalOutSampleR2(n) = results.allOutSampleR2AtAllThresholds(floor(results.allOptimalThresholds(n)/0.1), n);
        optimalOutSampleErrors(n) = results.allErrorsAtAllThresholds(floor(results.allOptimalThresholds(n)/0.1), n);
    end
    
    figure();
    plot(optimalOutSampleR2)
    hold on;
    plot(optimalOutSampleErrors)
elseif isnan(settings.Threshold)
    % Compare synthetic and cross-validated coefficients, and cumulative R^2 of
    % each taxon
    coefficientsComparison = table(results.syntheticCoefficients, results.crossValidatedCoefficients, results.cumulativeR2, 'VariableNames', {'syntheticCoefficients', 'crossValidatedCoefficients', 'cumulative R^2'});
end
end

% [~, idx] = sort(results.cumulativeR2, 'descend');
% selectedIdx = idx(1:30);
% recoveredCoefficients = zeros(100, 1);
% recoveredCoefficients(selectedIdx) = 1;
% [results.syntheticCoefficients, recoveredCoefficients]