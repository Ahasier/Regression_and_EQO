function [binaryCoefficients, optimalGroupSize] = runEQO(trainingData, trainingOutput, settings)
% Initialize variables
numTaxa = size(trainingData, 2);
% aicValues = zeros(1, numTaxa);
allAicValues = {[]};
allCoefficients = zeros(numTaxa, numTaxa);

% Save training data and output to CSV files for R to read
timestampStr = datestr(now, 'ddmmyy_HHMMSS');
trainingDataFilename = ['EQO/data/trainingData_', timestampStr, '.csv'];
trainingOutputFilename = ['EQO/data/trainingOutput_', timestampStr, '.csv'];

csvwrite(trainingDataFilename, trainingData);
csvwrite(trainingOutputFilename, trainingOutput);

% Loop across different group sizes
for maxGroupSize = 1:numTaxa
    % Call R script using system command
    coefficientsFilename = ['EQO/data/coefficients_', timestampStr, '_', num2str(maxGroupSize), '.csv'];
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
    allAicValues{groupSize}(end + 1) = computeAICofTheAssemblage(groupSize, trainingData, trainingOutput, coefficients, settings);
    
    % Store coefficients results
    allCoefficients(:, groupSize) = coefficients;
    
    % Delete the intermediate coefficients CSV file
    delete(coefficientsFilename);
end

% Get the average aic value at each actual group size
aicValues = averageAicValue(allAicValues);

% Find the optimal group size from minimum AIC value
[~, optimalGroupSize] = findMinimalAic(aicValues);

% Get the recovered coefficients corresponding to the optimal group size
binaryCoefficients = allCoefficients(:, optimalGroupSize);

% Delete the intermediate CSV files
delete(trainingDataFilename);
delete(trainingOutputFilename);
end

function aicValues = averageAicValue(allAicValues)
len = length(allAicValues);
aicValues = zeros(1, len);
for n = 1:len
    if ~isempty(allAicValues{n})
        aicValues(n) = NaN;
    else
        aicValues(n) = mean(allAicValues{n});
    end
end
end