function varargout = SetPathsForDataAndResults(varargin)
% Define the paths for data and results directories
dataPath = 'data/';
resultsPath = 'results/';
betaResultsPath = 'results/Betas/';
NAndXPath = [resultsPath,'NX/'];

% Initialize the output counter
outputCounter = 0;

% Loop through the input arguments
for i = 1:nargin
    if strcmp(varargin{i}, 'data')
        outputCounter = outputCounter + 1;
        varargout{outputCounter} = dataPath;
    elseif strcmp(varargin{i}, 'results')
        outputCounter = outputCounter + 1;
        varargout{outputCounter} = resultsPath;
    elseif strcmp(varargin{i}, 'betaResults')
        outputCounter = outputCounter + 1;
        varargout{outputCounter} = betaResultsPath;
    elseif strcmp(varargin{i}, 'NAndX')
        outputCounter = outputCounter + 1;
        varargout{outputCounter} = NAndXPath;
    else
        error('Unsupported input argument: %s', varargin{i});
    end
end
end
