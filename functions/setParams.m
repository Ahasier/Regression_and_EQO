function [numPermutations, phylogenyDependency, noiseLevel, meshGrid] = setParams(paramsFilename)
% Load configuration from JSON file
params = jsondecode(fileread(paramsFilename));

% Set parameters from configuration
numPermutations = params.numPermutations;
phylogenyDependency = params.phylogenyDependency;
noiseLevel = params.noiseLevel;
meshGrid = params.meshGrid;
end