clc; clear;

regressionMethod = 'LASSO';

% numberOfTaxaInAGroup_list = 5:5:50;
% numSamples_list = 10:10:200;

numberOfTaxaInAGroup_list = 10;
numSamples_list = 100;

% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations(regressionMethod);

% Create a full identifier string based on the options
fullIdentifier = createIdentifier(settings);

% Run computeAndSaveRegressionResults to get results files
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
function plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
% fullIdentifier2 = [fullIdentifier(1:end - 1),'ff'];
% fullIdentifier1 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_usePhylogenyOn';
% fullIdentifier2 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_usePhylogenyOff';
% fullIdentifier01 = [fullIdentifier,'_noise01'];
% fullIdentifier10 = [fullIdentifier,'_noise10'];
% 
% % load OLSMap and EQOMap from file
OLSMap = loadMapFromFile('OLS', fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMapNoise01 = loadMapFromFile('OLS', fullIdentifier01, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMapNoise10 = loadMapFromFile('OLS', fullIdentifier10, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% EQOMapNoise1 = loadMapFromFile('EQO', fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% EQOMapNoise01 = loadMapFromFile('EQO', fullIdentifier01, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% EQOMapNoise10 = loadMapFromFile('EQO', fullIdentifier10, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

% OLSMapUsePhylogeny = loadMapFromFile('OLS', fullIdentifier1, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% LASSOMapUsePhylogeny = loadMapFromFile('LASSO', fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% EQOMap = loadMapFromFile('EQO', fullIdentifier2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% NonlinearMap = loadMapFromFile('Nonlinear', fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

% fullIdentifier2 = [fullIdentifier(1:end - 2),'n'];

% % plot heatmaps
% figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
% subplot(2,3,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapNoise01, 'Noise level 0.1');
% subplot(2,3,2)
% h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapNoise1, 'Noise level 1');
% subplot(2,3,3)
% h3 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapNoise10, 'Noise level 10');
% subplot(2,3,4)
% h4 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise01, 'Noise level 0.1');
% subplot(2,3,5)
% h5 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise1, 'Noise level 1');
% subplot(2,3,6)
% h6 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise10, 'Noise level 10');

% subplot(2,4,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
% subplot(2,4,2)
% h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS');
% subplot(2,4,3)
% h3 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapUsePhylogeny, 'Phylo-regularized OLS');
% subplot(2,4,4)
% h4 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapUsePhylogeny, 'Phylo-regularized LASSO');
% subplot(2,4,5)
% h5 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapUsePhylogeny, 'Phylo-regularized OLS', OLSMap, 'OLS');
% subplot(2,4,6)
% h6 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapUsePhylogeny, 'Phylo-regularized LASSO', OLSMap, 'OLS');
% subplot(2,4,7)
% h6 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMapUsePhylogeny, 'Phylo-regularized OLS', EQOMap, 'EQO');
% subplot(2,4,8)
% h6 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapUsePhylogeny, 'Phylo-regularized LASSO', EQOMap, 'EQO');
% saveas(gcf,['results/plots/Heatmap',fullIdentifier,'_noiseCompare.jpg'])

figure();
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS on Phylo-independent data',1);
end

function Map = loadMapFromFile(regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
Map = accessToAccMap('results/Accuracy/', regressionMethod, fullIdentifier, 'Load');
Map = Map(numberOfTaxaInAGroup_list/meshGrid.TaxaGroup, numSamples_list/meshGrid.Samples);
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