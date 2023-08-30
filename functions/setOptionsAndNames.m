function [settings, fullIdentifier] = setOptionsAndNames(varargin)
% SETOPTIONSANDNAES configures options and creates a unique identifier string.
% INPUTS:
%   varargin: Name-value pairs for setting various options
% OUTPUTS:
%   settings: A structure containing the configuration settings
%   fullIdentifier: A unique identifier string based on the options

% Confirm that an even number of input arguments have been provided.
% Otherwise, it’s an error.
if mod(length(varargin), 2) ~= 0
    error('Expected input arguments in name-value pairs.');
end

% Initialize options as an empty struct
settings = struct();

% Populate the settings struct with name-value pairs
for i = 1:2:length(varargin)
    % varargin{i} is the option, varargin{i+1} is the value
    settings.(varargin{i}) = varargin{i+1};
end

% Initialize fullIdentifier as an empty string
fullIdentifier = [];

% Create a full identifier string based on the options
optionnames = fieldnames(settings);
for n = 1:length(optionnames)
    strName = getShortName(optionnames{n});
    strValue = getValueAsString(settings.(optionnames{n}));
    addStr = ['_', strName, strValue];
    fullIdentifier = [fullIdentifier, addStr];
end

% Set default value for BetaEps if it's not provided
if ~isfield(settings, 'BetaEps')
    settings.BetaEps = 0;
end
end

%% Helper functions
function strName = getShortName(optionName)
% Convert full option names to shorter versions or an empty string
if strcmp(optionName, 'Threshold')
    strName = 'Tr';
elseif strcmp(optionName, 'L0Option')
    strName = [];
else
    strName = optionName;
end
end

function strValue = getValueAsString(value)
% Convert the value to a string, regardless of its original type
if isstring(value)
    strValue = value;
else
    strValue = num2str(value);
end
end