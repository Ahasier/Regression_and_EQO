function existingData = loadOrInitializeAccuracyResultsFile(filename)
% Check if the results file already exists
if isfile(filename)
    % If file exists, read the existing data from the CSV file
    existingData = csvread(filename);
else
    % If file does not exist, initialize a new dataset with zeros
    existingData = 0;
    % Write the initialized dataset to a new CSV file
    csvwrite(filename, existingData);
end
end