function readAndPlotTCM(numberOfTaxaInAGroup_list, numSamples_list, regressionMethod, fullIdentifier)
% READANDPLOTTCM reads and plots the TCM based on the provided parameters.
% Input:
% - numberOfTaxaInAGroup_list: List of numbers indicating the number of taxa in each group.
% - numSamples_list: List of sample sizes.
% - regressionMethod: The regression method used.
% - fullIdentifier: An identifier for the dataset or method.

% Define the indices for taxa and assemblages
taxaIndices = 1:100;
assemblageIndices = 1:100;

% Create a new figure with specified renderer and position
figure('Renderer', 'painters', 'Position', [0 0 1000 1000]);

% Loop through each taxa group size
for m = 1:length(numSamples_list)
    numSamples = numSamples_list(m);
    
    % Loop through each sample size
    for n = 1:length(numberOfTaxaInAGroup_list)
        numberOfTaxaInAGroup = numberOfTaxaInAGroup_list(n);
        
        % Read the TCM table for the given parameters
        [TCM, R2] = readTCMTable(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples);
        
        % Read the ground truth from results for the given parameters
        results = loadResults(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples);
        groundTruth = results.syntheticCoefficients > 0;
        binaryResults = results.crossValidatedCoefficients;
        importanceValue = results.importanceValues;
        
        % Arrange TCM and ground truth in descending order
        [arrangedTCM, arrangedeImportanceValue, arrangedBinaryResults, arrangedGroundTruth, arrangedR2] = arrangeTCM(TCM, binaryResults, importanceValue, groundTruth, R2);
        
        % Plot the TCM and R^2 values
        plotTCM(arrangedTCM, importanceValue, binaryResults, arrangedR2, groundTruth, regressionMethod, taxaIndices, assemblageIndices, n, m, numberOfTaxaInAGroup, numSamples, numberOfTaxaInAGroup_list, numSamples_list);
    end
end

% Save the figure as a .jpg file in the specified directory
saveas(gcf, ['results/plots/realdata_standardized/TCM', regressionMethod, fullIdentifier, '.jpg'])
end

%% Helper functions
function [TCM, R2] = readTCMTable(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples)
global paths

% Read the TCM table from saved CSV file
tcmFilename = [paths.betaResultsPath, 'realdata_standardized/Betas_', regressionMethod, fullIdentifier, '_nSpl', num2str(numSamples), '.mat'];
% TCMtable = csvread(tcmFilename, 1, 0);
load(tcmFilename, 'results');
TCMtable = results.TCM;

% Extract the TCM values (all rows except the last one)
TCM = TCMtable{1:end - 1, :};

% Extract the R^2 values (last row)
R2 = TCMtable{end, :};
end

function [arrangedTCM, arrangedeImportanceValue, arrangedBinaryResults, arrangedGroundTruth, arrangedR2] = arrangeTCM(TCM, binaryResults, importanceValue, groundTruth, R2)
[~, arrangeOrder] = sort(groundTruth, 'descend');
arrangedGroundTruth = groundTruth(arrangeOrder, :);
arrangedeImportanceValue = importanceValue(arrangeOrder, :);
arrangedBinaryResults = binaryResults(arrangeOrder, :);

[~, OrderOfR2] = sort(R2,'descend');
arrangedTCM = TCM(:, OrderOfR2);
arrangedR2 = R2(OrderOfR2);
end

function plotTCM(TCM, importanceValue, binaryResults, R2, groundTruth, regressionMethod, taxaIndices, assemblageIndices, n, m, numberOfTaxaInAGroup, numSamples, numberOfTaxaInAGroup_list, numSamples_list)
% Define the number of rows and columns for your subplots
nh = length(numSamples_list); % number of rows
nw = length(numberOfTaxaInAGroup_list); % number of columns

% Define margins for spacing around subplots
marginWidth = 0.06;
marginHeight = 0.04;
gapBetweenPlots = 0.04; % Additional gap between R^2 and TCM plots
marginWidthGroundTruth = 0.01;

% Adjusted width and height for TCM, R2, and groundTruth subplots
widthGroundTruth = 0.01; % Narrow width for groundTruth heatmap
widthTCM = (1 - (nw + 4) * marginWidth - 2 * widthGroundTruth) / nw; % Adjusted width for TCM
heightTCM = (1 - (1.4 * nh) * marginHeight - gapBetweenPlots * (nh-1)) / (1.4 * nh);
heightR2 = heightTCM / 20;  % Make R2 height 1/20th of TCM height

% Calculate the left and bottom positions for the subplots
leftTCM = (n-1) * (widthTCM + 3 * widthGroundTruth + 5 * marginWidthGroundTruth + marginWidth) + marginWidth; % Position for importanceValue plot
leftImportanceValue = leftTCM + widthTCM + marginWidthGroundTruth; % Position for TCM plot
leftBinaryResults = leftImportanceValue + widthGroundTruth + marginWidthGroundTruth; % Position for binaryResults plot
leftGroundTruth = leftBinaryResults + widthGroundTruth + marginWidthGroundTruth; % Position for groundTruth plot
bottomR2 = 1 - (m * (heightTCM + heightR2 + 2 * marginHeight + gapBetweenPlots));
bottomTCM = bottomR2 + heightR2 + marginHeight;

% Create the TCM subplot
axTCM = axes('Position', [leftTCM, bottomTCM, widthTCM, heightTCM]);
h = imagesc(assemblageIndices, taxaIndices, TCM);
title(axTCM, ['TCM (K' num2str(numberOfTaxaInAGroup) ', nSpl' num2str(numSamples) ')']);
if n == 1
    ylabel(axTCM, "Taxa Index");
end

% Set the missing data pixels to grey
h.AlphaData = ones(size(TCM));
h.AlphaData([7 14 19 23 25 26 30 33 42 46 49 59 61 70 73 80 82 88], :) = 0.2;
h.AlphaData([12 16 18 32 34 54 66 84 94 96 99], :) = 0.5;
h.AlphaData([6 13 15 24 27 35 36 44 50 55 78 85 91], :) = 0.7;
ax = gca; ax.Color = [1 1 1];

% Create the importanceValue subplot
axImportanceValue = axes('Position', [leftImportanceValue, bottomTCM, widthGroundTruth, heightTCM]);
imagesc(1, taxaIndices, repmat(importanceValue, 1, size(TCM, 2)));
title(axImportanceValue, '$\eta$', 'Interpreter', 'latex');
set(axImportanceValue, 'XTick', [], 'YTick', []); % Hide ticks

% Create the binaryResults subplot
axBinaryResults = axes('Position', [leftBinaryResults, bottomTCM, widthGroundTruth, heightTCM]);
imagesc(1, taxaIndices, repmat(binaryResults, 1, size(TCM, 2)));
title(axBinaryResults, '$\hat{\beta}$', 'Interpreter', 'latex');
set(axBinaryResults, 'XTick', [], 'YTick', []); % Hide ticks

% Create the groundTruth subplot
axGroundTruth = axes('Position', [leftGroundTruth, bottomTCM, widthGroundTruth, heightTCM]);
imagesc(1, taxaIndices, groundTruth);
title(axGroundTruth, '$\beta_0$', 'Interpreter', 'latex');
set(axGroundTruth, 'XTick', [], 'YTick', []); % Hide ticks

% Create the R2 subplot
axR2 = axes('Position', [leftTCM, bottomR2, widthTCM, heightR2]);
imagesc(assemblageIndices, 1, R2);
caxis([0 1]);
if n == 1
ylabel(axR2, 'R^2');
end
if m == nh
xlabel(axR2, 'Assemblage Index');
end
set(axR2, 'XTick', [], 'YTick', []); % Hide ticks

% Apply a consistent style across all plots
% Here you can define a function 'plotstyle' to apply styles like fonts, colors, etc.
plotstyle(axImportanceValue, 1);
plotstyle(axTCM, 1);
plotstyle(axBinaryResults, 1);
plotstyle(axGroundTruth, 1);
plotstyle(axR2, 1);

end