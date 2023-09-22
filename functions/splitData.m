function [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput, varargin)
totalSamples = size(abundanceData, 1);
idx = randperm(totalSamples);

% Assuming a 50-50 split for training and testing
splitPoint = floor(0.5 * totalSamples);

trainingIdx = idx(1:splitPoint);
testIdx = idx(splitPoint+1:end);

trainingData = abundanceData(trainingIdx, :);
testData = abundanceData(testIdx, :);

trainingOutput = functionalOutput(trainingIdx);
testOutput = functionalOutput(testIdx);
end