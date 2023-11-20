numberOfTaxaInAGroup_list = 5:5:50;
numSamples_list = 10:10:200;

% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations('OLS');

% % Figure 2: Comparison between OLS and EQO on binary ground truth
fullIdentifierF2 = '_Beta01_BetaEps0_TrNaN';
EQOMap = loadMapFromFile('EQO', fullIdentifierF2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
OLSMap = loadMapFromFile('OLS', fullIdentifierF2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
subplot(1,3,1)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS');
subplot(1,3,2)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
subplot(1,3,3)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS', EQOMap, 'EQO');
% 
% % Figure 3: Accuracy Change on Non-binary Ground Truth
% fullIdentifierF3 = '_Beta01_BetaEps0.5_TrNaN';
% EQOMap = loadMapFromFile('EQO', fullIdentifierF3, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMap = loadMapFromFile('OLS', fullIdentifierF3, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
% subplot(1,3,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS');
% subplot(1,3,2)
% h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
% subplot(1,3,3)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS', EQOMap, 'EQO');

% Figure 4: Comparison between OLS and EQO on non-binary ground truth and
% mimic-real abundance data
% fullIdentifierF4 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn';
% EQOMap = loadMapFromFile('EQO', fullIdentifierF4, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMap = loadMapFromFile('OLS', fullIdentifierF4, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
% subplot(1,3,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS');
% subplot(1,3,2)
% h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
% subplot(1,3,3)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap, 'OLS', EQOMap, 'EQO');

% Figure 5: Comparison between LASSO, OLS and EQO
% fullIdentifierF5 = '_Beta01_BetaEps0_TrNaN';
% fullIdentifierF5_2 = '_Beta01_BetaEps0_TrNaN_maxLambda10';
% EQOMap = loadMapFromFile('EQO', fullIdentifierF5, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMap = loadMapFromFile('OLS', fullIdentifierF5, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% LASSOMap = loadMapFromFile('LASSO', fullIdentifierF5_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% figure('Renderer', 'painters', 'Position', [0 0 1500 320]);
% subplot(1,3,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap, 'LASSO');
% subplot(1,3,2)
% h2 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap, 'LASSO', OLSMap, 'OLS');
% subplot(1,3,3)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap, 'LASSO', EQOMap, 'EQO');

% Figure 6: Comparison between phylo-regularized OLS and EQO on 
% phylogenetcally-dependent non-binary ground truth and mimic-real 
% abundance data
% fullIdentifierF6_1 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_usePhylogenyOn';
% fullIdentifierF6_2 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_usePhylogenyOff';
% fullIdentifierF6_3 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_maxLambda10_usePhylogenyOn';
% EQOMap2 = loadMapFromFile('EQO', fullIdentifierF6_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMap2 = loadMapFromFile('OLS', fullIdentifierF6_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% OLSMap1 = loadMapFromFile('OLS', fullIdentifierF6_1, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% LASSOMap1 = loadMapFromFile('LASSO', fullIdentifierF6_3, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% figure('Renderer', 'painters', 'Position', [0 0 1500 500]);
% subplot(2,4,1)
% h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap2, 'OLS');
% subplot(2,4,2)
% h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap1, 'phylo-regularized OLS');
% subplot(2,4,3)
% h3 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap1, 'phylo-regularized LASSO');
% subplot(2,4,4)
% h3 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap2, 'EQO');
% subplot(2,4,5)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap1, 'phylo-regularized OLS', OLSMap2, 'OLS');
% subplot(2,4,6)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap1, 'phylo-regularized LASSO', OLSMap2, 'OLS');
% subplot(2,4,7)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, OLSMap1, 'phylo-regularized OLS', EQOMap2, 'EQO');
% subplot(2,4,8)
% h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMap1, 'phylo-regularized LASSO', EQOMap2, 'EQO');


%% Helper functions
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