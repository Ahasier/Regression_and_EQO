%% Add path
added_path = [pwd,'/functions'];
addpath(added_path);

%% Set parameters and input arguments
numPermutations = 10; % Average over 10 permutations
phylogenyDependency = 0; % Groung truth beta's are Phylogenetically irrelavent
noiseLevel = 1; % Noise level = 1

numberOfTaxaInAGroup_list = [10 30 50]; % Set the number of taxa within a syntetic functional group
numSamples_list = [40 80 120]; % Set the number of samples

meshGrid.TaxaGroup = 5; % The results map use a mesh grid of 5:5:50 numberOfTaxaInAGroup
meshGrid.Samples = 10; % The results map use a mesh grid of 10:10:200 numSamples

regressionMethod = 'OLS';
varargin = {'Beta0', 1, 'BetaEps', 0.5, 'Threshold', 'cv','RealAbd','On'};
[settings, fullIdentifier] = setOptionsAndNames(varargin{:});

%% Run computeAndSaveRegressionResults to get results files
for numberOfTaxaInAGroup = numberOfTaxaInAGroup_list
    for numSamples = numSamples_list
        computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, varargin{:})
    end
end

%% Plot heatmaps
plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

%% Helper functions
function plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
% load OLSBeta0Map and EQOMap from file
OLSBeta0Map = accessToAccMap(regressionMethod, fullIdentifier, 'Load');
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

saveas(gcf,['results/Heatmap',fullIdentifier,'.jpg'])
end

function resultMap = accessToAccMap(regressionMethod, fullIdentifier, saveOrLoad)
filename = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end

function resultMap = accessToEQOMap(Beta0, BetaEps, saveOrLoad)
filename = ['results/AccMCC_Beta0',num2str(Beta0),'_BetaEps',num2str(BetaEps),'_RealAbdOn_EQO.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end