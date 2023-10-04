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
for m = 1:length(numberOfTaxaInAGroup_list)
    numberOfTaxaInAGroup = numberOfTaxaInAGroup_list(m);
    
    % Loop through each sample size
    for n = 1:length(numSamples_list)
        numSamples = numSamples_list(n);
        
        % Read the TCM table for the given parameters
        [TCM, R2] = readTCMTable(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples);
        
        % Read the ground truth from results for the given parameters
        results = loadResults(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples);
        groundTruth = results.syntheticCoefficients > 0;
        
        % Plot the TCM and R^2 values
        plotTCM(TCM, R2, groundTruth, regressionMethod, taxaIndices, assemblageIndices, n, m, numberOfTaxaInAGroup, numSamples, numberOfTaxaInAGroup_list, numSamples_list);
    end
end

% Save the figure as a .jpg file in the specified directory
saveas(gcf, ['results/plots/TCM', regressionMethod, fullIdentifier, '.jpg'])
end

%% Helper functions
function [TCM, R2] = readTCMTable(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples)
global paths

% Read the TCM table from saved CSV file
tcmFilename = [paths.tcmResultsPath, 'TCM_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.csv'];
TCMtable = csvread(tcmFilename, 1, 0);

% Extract the TCM values (all rows except the last one)
TCM = TCMtable(1:end - 1, :);

% Extract the R^2 values (last row)
R2 = TCMtable(end, :);
end

function plotTCM(TCM, R2, groundTruth, regressionMethod, taxaIndices, assemblageIndices, n, m, numberOfTaxaInAGroup, numSamples, numberOfTaxaInAGroup_list, numSamples_list)
% Define the number of rows and columns for your subplots
nh = length(numberOfTaxaInAGroup_list); % number of rows
nw = length(numSamples_list); % number of columns

% Define margins for spacing around subplots
marginWidth = 0.1;
marginHeight = 0.04;
gapBetweenPlots = 0.04; % Additional gap between R^2 and TCM plots
marginWidthGroundTruth = 0.02;

% Adjusted width and height for TCM, R2, and groundTruth subplots
widthGroundTruth = 0.01; % Narrow width for groundTruth heatmap
widthTCM = (1 - (nw + 1) * marginWidth - widthGroundTruth) / nw; % Adjusted width for TCM considering groundTruth width
heightTCM = (1 - (1.4 * nh) * marginHeight - gapBetweenPlots * (nh-1)) / (1.4 * nh);
heightR2 = heightTCM / 20;  % Make R2 height 1/20th of TCM height

% Calculate the left and bottom positions for the subplots
leftTCM = (n-1) * (widthTCM + marginWidth + widthGroundTruth) + marginWidth; % Position for TCM plot
leftGroundTruth = leftTCM + widthTCM + marginWidthGroundTruth; % Position for groundTruth plot to the right of TCM
bottomR2 = 1 - (m * (heightTCM + heightR2 + 2 * marginHeight + gapBetweenPlots));
bottomTCM = bottomR2 + heightR2 + marginHeight;

% Create the TCM subplot
ax1 = axes('Position', [leftTCM, bottomTCM, widthTCM, heightTCM]);
imagesc(assemblageIndices, taxaIndices, TCM);

% Adjust the title position for TCM
titleObjTCM = title(['TCM (', regressionMethod, ', K', num2str(numberOfTaxaInAGroup), ', nSpl', num2str(numSamples), ')'],  'fontsize', 18);
titlePosTCM = get(titleObjTCM, 'Position');
set(titleObjTCM, 'Position', [titlePosTCM(1) - 6, titlePosTCM(2), titlePosTCM(3)]); % Adjust the value '5' as needed

if n == 1
    ylabel("Taxa index", 'fontsize', 18);
end
ax1 = plotstyle(gca);

% Create the groundTruth subplot
ax0 = axes('Position', [leftGroundTruth, bottomTCM, widthGroundTruth, heightTCM]);
imagesc(1, taxaIndices, groundTruth);

% Hide the ticks for groundTruth
set(ax0, 'XTick', [], 'YTick', []);

% Adjust the title position for groundTruth
titleObj = title('Ground Truth', 'fontsize', 18);
titlePos = get(titleObj, 'Position');
set(titleObj, 'Position', [titlePos(1) + 1.5, titlePos(2), titlePos(3)]);

ax0 = plotstyle(gca);

% Create the R2 subplot
ax2 = axes('Position', [leftTCM, bottomR2, widthTCM, heightR2]);
imagesc(assemblageIndices, 1, R2);
if n == 1
    ylabel("R^2", 'fontsize', 18);
end
if m == 3
    xlabel("assemblage index", 'fontsize', 18);
end
ax2 = plotstyle(gca);
end