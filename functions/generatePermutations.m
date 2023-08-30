function [trainingSet, testSet] = generatePermutations(numberOfSamples)
trainingSet = randperm(numberOfSamples, numberOfSamples / 2);
testSet = 1:numberOfSamples;
testSet(trainingSet) = [];
end