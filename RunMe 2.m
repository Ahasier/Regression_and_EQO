clc;clear
%% Define global variable paths for where to load data or store results
global paths

%% Add necessary paths
initializePaths();

% Add the Path to `Rscript` to PATH in MATLAB
setenv('PATH', [getenv('PATH') ':/usr/local/bin/']);

%% Set parameters and input arguments
numPermutations = 1; % Average over 100 permutations
phylogenyDependency = 0; % Groung truth beta's are Phylogenetically irrelavent
noiseLevel = 1; % Noise level = 1

numberOfTaxaInAGroup_list = [10 30 50]; % Set the number of taxa within a syntetic functional group
numSamples_list = [40 80 120]; % Set the number of samples

meshGrid.TaxaGroup = 5; % The results map use a mesh grid of 5:5:50 numberOfTaxaInAGroup
meshGrid.Samples = 10; % The results map use a mesh grid of 10:10:200 numSamples

regressionMethod = 'EQO';
varargin = {'Beta0', 1, 'BetaEps', 0.5, 'Threshold', 0.5, 'RealAbd','On', 'requirePositivity', 'On'};

[settings, fullIdentifier] = setOptionsAndNames(varargin{:});
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

%% Run computeAndSaveRegressionResults to get results files
for numberOfTaxaInAGroup = numberOfTaxaInAGroup_list
    for numSamples = numSamples_list
        results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier);
    end
end

%% Plot heatmaps
plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

%% Helper functions
function initializePaths()
functionsPath = [pwd, '/functions'];
addpath(functionsPath);
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