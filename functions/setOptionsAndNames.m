function [settings, fullIdentifier] = setOptionsAndNames()
% SETOPTIONSANDNAES configures options and creates a unique identifier string.
% INPUTS:
%   varargin: Name-value pairs for setting various options
% OUTPUTS:
%   settings: A structure containing the configuration settings
%   fullIdentifier: A unique identifier string based on the options

% Load the configuration settings from config.json
fid = fopen('configurations/config.json', 'r');
rawData = fread(fid, inf, '*char')';
fclose(fid);
settings = jsondecode(rawData);

% Convert string "NaN" to MATLAB NaN
if isfield(settings, 'Threshold') && strcmp(settings.Threshold, 'NaN')
    settings.Threshold = NaN;
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