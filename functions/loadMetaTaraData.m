% Load Meta Tara data
function metaTaraData = loadMetaTaraData(pathToData)
pathToMetaTaraData = [pathToData, 'Meta_Tara.csv'];
metaTaraData = readtable(pathToMetaTaraData);
end