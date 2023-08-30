function fullPath = getCompleteAccuracyFilePath(regressionMethod, fullIdentifier, resultsPath)
fullPath = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];
end