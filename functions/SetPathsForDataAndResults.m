function paths = SetPathsForDataAndResults(varargin)
% Define the paths for data and results directories
dataPath = 'data/';
resultsPath = 'results/';
betaResultsPath = 'results/Betas/';
NAndXPath = [resultsPath,'NX/'];

% Loop through the input arguments
for i = 1:nargin
    if strcmp(varargin{i}, 'data')
        paths.data = dataPath;
    elseif strcmp(varargin{i}, 'results')
        paths.resultsPath = resultsPath;
    elseif strcmp(varargin{i}, 'betaResults')
        paths.betaResultsPath = betaResultsPath;
    elseif strcmp(varargin{i}, 'NAndX')
        paths.NAndX = NAndXPath;
    else
        error('Unsupported input argument: %s', varargin{i});
    end
end
end
