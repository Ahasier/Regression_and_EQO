function changeRecoveredBetaUsingCappedImportanceValues(regressionMethod, betaEps, realAbd, usePhylogeny)
% Initialize all neccessary parametters from configurations files
[~, ~, ~, ~, settings] = initializations(regressionMethod);

settings.BetaEps = betaEps;

if ~strcmp(realAbd, 'EMPTY_ARRAY')
    settings.RealAbd = realAbd;
end

if ~strcmp(usePhylogeny, 'EMPTY_ARRAY')
    settings.usePhylogeny = usePhylogeny;
end

% Create a full identifier string based on the options
fullIdentifier = createIdentifier(settings);

numberOfTaxaInAGroup_list = 5:5:50;
numSamples_list = 10:10:200;

global paths

if ~usingRealData(settings)
    for n = 1:length(numberOfTaxaInAGroup_list)
        numberOfTaxaInAGroup = numberOfTaxaInAGroup_list(n);
        for m = 1:length(numSamples_list)
            numSamples = numSamples_list(m);
            betaResultsFilePath = [paths.betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];
            if isfile(betaResultsFilePath)
                load(betaResultsFilePath, 'results');
                
                results.crossValidatedCoefficients = results.importanceValues > 0.5;
                
                results.accuracy = calculateAccuracy(results.crossValidatedCoefficients, results.syntheticCoefficients, settings);
                
                % Construct the file path where other simulation results will be saved
                betaResultsFilePath = [paths.betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];
                
                % Save various results to a MAT file for later analysis
                save(betaResultsFilePath, 'results');
            end
        end
    end
end
end