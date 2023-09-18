% Load tree data from a file
function treeData = loadRealTreeData()
global paths

filePath = [paths.data, 'tree100taxaReal.mat'];
treeData = load(filePath);
end