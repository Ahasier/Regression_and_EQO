clc;clear
%% Set parameters and input arguments
[numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup_list, numSamples_list, meshGrid, regressionMethod] = setParams();
[settings, fullIdentifier] = setOptionsAndNames();

%% Set paths
% Add necessary paths
initializePaths();

% Define global variable paths for where to load data or store results
global paths
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

% If using EQO, add the path to `Rscript` to PATH in MATLAB
if strcmp(regressionMethod, 'EQO')
    setenv('PATH', [getenv('PATH') ':/usr/local/bin/']);
end

%% Run computeAndSaveRegressionResults to get results files
for numberOfTaxaInAGroup = numberOfTaxaInAGroup_list
    for numSamples = numSamples_list
        results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier);
    end
end

%% Plot TCMs
readAndPlotTCM(numberOfTaxaInAGroup_list, numSamples_list, regressionMethod, fullIdentifier)

%% Plot heatmaps
plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

%% Helper functions
function initializePaths()
functionsPath = [pwd, '/functions/'];
addpath(functionsPath);
end

function [numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup_list, numSamples_list, meshGrid, regressionMethod] = setParams()
numPermutations = 100; % Average over 100 permutations
phylogenyDependency = 0; % Groung truth beta's are Phylogenetically irrelavent
noiseLevel = 1; % Noise level = 1

numberOfTaxaInAGroup_list = [10 30 50]; % Set the number of taxa within a syntetic functional group
numSamples_list = [40 80 100]; % Set the number of samples

meshGrid.TaxaGroup = 5; % The results map use a mesh grid of 5:5:50 numberOfTaxaInAGroup
meshGrid.Samples = 10; % The results map use a mesh grid of 10:10:200 numSamples

regressionMethod = 'OLS';
end

function plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
% load OLSBeta0Map and EQOMap from file
OLSBeta0Map = accessToAccMap('results/Accuracy/', regressionMethod, fullIdentifier, 'Load');
OLSBeta0Map = OLSBeta0Map(numberOfTaxaInAGroup_list/meshGrid.TaxaGroup, numSamples_list/meshGrid.Samples);

EQOMap = accessToEQOMap(settings.Beta0, settings.BetaEps, 'Load');
EQOMap = EQOMap(numberOfTaxaInAGroup_list/meshGrid.TaxaGroup, numSamples_list/meshGrid.Samples);

% plot heatmaps
figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
subplot(1,3,1)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
subplot(1,3,2)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSBeta0Map, 'OLS');
subplot(1,3,3)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSBeta0Map, 'OLS', EQOMap, 'EQO');

saveas(gcf,['results/plots/Heatmap',fullIdentifier,'.jpg'])
end

function resultMap = accessToAccMap(resultsPath, regressionMethod, fullIdentifier, saveOrLoad)
filename = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end

function resultMap = accessToEQOMap(Beta0, BetaEps, saveOrLoad)
filename = ['results/Accuracy/AccMCC_Beta0',num2str(Beta0),'_BetaEps',num2str(BetaEps),'_RealAbdOn_EQO.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end