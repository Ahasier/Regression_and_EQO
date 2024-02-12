function [binaryCoefficients, aicValue, optimalGroupSize] = runEQO(trainingData, trainingOutput, numberOfTaxaInAGroup, settings)
% Initialize variables
numTaxa = size(trainingData, 2);
% aicValues = zeros(1, numTaxa);
allAicValues = arrayfun(@(x) [], cell(numTaxa,1), 'UniformOutput', false);
allCoefficients = arrayfun(@(x) [], cell(numTaxa,1), 'UniformOutput', false);
lenAIC = min(numTaxa, length(trainingOutput) - 2);

% maxGroupSize = numberOfTaxaInAGroup + 15;
numSamples = 2 * length(trainingOutput);

if isfield(settings, "usePhylogeny")
    usePhylogenyIdentifier = settings.usePhylogeny;
else
    usePhylogenyIdentifier = [];
end

if isfield(settings, "RealAbd")
    RealAbdIdentifier = settings.RealAbd;
else
    RealAbdIdentifier = [];
end

% Save training data and output to CSV files for R to read
timestampStr = datestr(now, 'ddmmyy_HHMMSS');
trainingDataFilename = ['EQO/data/trainingData_', num2str(settings.BetaEps), '_', RealAbdIdentifier, '_',  usePhylogenyIdentifier, '_', num2str(numberOfTaxaInAGroup), '_',num2str(numSamples), '_', timestampStr, '.csv'];
trainingOutputFilename = ['EQO/data/trainingOutput_', num2str(settings.BetaEps), '_', RealAbdIdentifier, '_',  usePhylogenyIdentifier, '_', num2str(numberOfTaxaInAGroup), '_',num2str(numSamples), '_', timestampStr, '.csv'];

csvwrite(trainingDataFilename, trainingData);
csvwrite(trainingOutputFilename, trainingOutput);

% Loop across different group sizes
for maxGroupSize = 1:numberOfTaxaInAGroup + 20
    
    % Call R script using system command
    coefficientsFilename = ['EQO/data/coefficients_', num2str(settings.BetaEps), '_', RealAbdIdentifier, '_',  usePhylogenyIdentifier, '_', num2str(numberOfTaxaInAGroup), '_',num2str(numSamples), '_', timestampStr, '_', num2str(maxGroupSize), '.csv'];
    commandStr = sprintf('Rscript EQO/EQO_GA_script.R %d %s %s %s', maxGroupSize, trainingDataFilename, trainingOutputFilename, coefficientsFilename);
    system(commandStr);
    
    % Check for the existence of the coefficients.csv file
    maxAttempts = 120; % Maximum number of attempts to check for the file
    delayInSeconds = 10; % Delay between checks
    
    for attempt = 1:maxAttempts
        if exist(coefficientsFilename, 'file') == 2
            % File exists, break out of the loop
            break;
        else
            % If file doesn't exist, pause and then check again
            pause(delayInSeconds);
        end
    end
    
    % Read results back from R
    coefficients = csvread(coefficientsFilename, 1, 0);
    
    % Count the actual group size
    groupSize = sum(coefficients > 0);
    
    % Calculate the AIC value in this group size.
    % aicValue = computeAIC(groupSize, trainingData, trainingOutput, coefficients);
    
    if groupSize < lenAIC
        % Calculate the AIC value in this group size.
        uncorrectedAic = computeAIC(groupSize, trainingData, trainingOutput, coefficients);
        allAicValues{groupSize}(end + 1) = correctAIC(uncorrectedAic, groupSize, numSamples);
    end
    
    % Store coefficients results
    allCoefficients{groupSize}(:,end + 1) = coefficients;
    
    % Delete the intermediate coefficients CSV file
    delete(coefficientsFilename);
end

% Get the average aic value at each actual group size
aicValues = averageAicValue(allAicValues);

% Find the optimal group size from minimum AIC value
[aicValue, optimalGroupSize] = findMinimalAic(aicValues);

% Get the recovered coefficients corresponding to the optimal group size
binaryCoefficients = getOptimalCoefficients(allCoefficients, optimalGroupSize);

% Delete the intermediate CSV files
delete(trainingDataFilename);
delete(trainingOutputFilename);
end

function aicValues = averageAicValue(allAicValues)
len = length(allAicValues);
aicValues = zeros(1, len);
for n = 1:len
    if isempty(allAicValues{n})
        aicValues(n) = NaN;
    else
        aicValues(n) = mean(allAicValues{n});
    end
end
end

function binaryCoefficients = getOptimalCoefficients(allCoefficients, optimalGroupSize)
coefficientsOfOptimalGroupSize = allCoefficients{optimalGroupSize};

freqOfTaxa = sum(coefficientsOfOptimalGroupSize, 2);
[~, idx] = sort(freqOfTaxa, 'descend');

numTaxa = size(coefficientsOfOptimalGroupSize, 1);

binaryCoefficients = zeros(numTaxa, 1);
binaryCoefficients(idx(1:optimalGroupSize)) = 1;
end