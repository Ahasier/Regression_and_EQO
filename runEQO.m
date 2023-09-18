function [binaryCoefficients, optimalGroupSize] = runEQO(trainingData, trainingOutput)
% Initialze variables
numTaxa = size(trainingData, 2);
aicValues = zeros(1, numTaxa);
allCoefficients = zeros(numTaxa, numTaxa);

% Save training data and output to CSV files for R to read
csvwrite('EQO/data/trainingData.csv', trainingData);
csvwrite('EQO/data/trainingOutput.csv', trainingOutput);

% Loop across different group sizes
for groupSize = 1:numTaxa
    % Add the Path to `Rscript` to PATH in MATLAB
    setenv('PATH', [getenv('PATH') ':/usr/local/bin/']);
    
    % Call R script using system command
    commandStr = sprintf('Rscript EQO/EQO_GA_script.R %d', groupSize);
    system(commandStr);
    
    % Check for the existence of the coefficients.csv file
    maxAttempts = 120; % Maximum number of attempts to check for the file
    delayInSeconds = 10; % Delay between checks
    
    for attempt = 1:maxAttempts
        if exist('data/coefficients.csv', 'file') == 2
            % File exists, break out of the loop
            break;
        else
            % If file doesn't exist, pause and then check again
            pause(delayInSeconds);
        end
    end
    
    % Read results back from R
    coefficients = csvread('data/coefficients.csv', 1, 0);
    
    % Sort coefficients by descend for later use
    [~, sortedTaxaIndices] = sort(coefficients, 'descend');
    
    % calculate the AIC value in this group size.
    aicValues(groupSize) = computeAIC(groupSize, numTaxa, trainingData, trainingOutput, sortedTaxaIndices);
    
    % Store coefficients results
    allCoefficients(:, groupSize) = coefficients;
end

% Find the optimal group size from minimum AIC value
[~, optimalGroupSize] = findMinimalAic(aicValues);

% Get the recovered coefficients corresponding to the optimal group size
binaryCoefficients = allCoefficients(:, optimalGroupSize);

% Delete the intermediate CSV files
delete('EQO/data/trainingData.csv');
delete('EQO/data/trainingOutput.csv');
delete('data/coefficients.csv');
end
