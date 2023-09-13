% Load tree data from a file
function treeData = loadRealTreeData()
% Set path to data
pathToData = SetPathsForDataAndResults('data');
% Load tree from file
filePath = [pathToData, 'tree100taxaReal.mat'];
treeData = load(filePath);
end