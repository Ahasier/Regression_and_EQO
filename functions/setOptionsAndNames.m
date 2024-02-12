function settings = setOptionsAndNames()
% SETOPTIONSANDNAES configures options and creates a unique identifier string.
% INPUTS:
%   varargin: Name-value pairs for setting various options
% OUTPUTS:
%   settings: A structure containing the configuration settings
%   fullIdentifier: A unique identifier string based on the options

% Load the configuration settings from config.json
fid = fopen('configurations/settingsConfig.json', 'r');
rawData = fread(fid, inf, '*char')';
fclose(fid);
settings = jsondecode(rawData);

% Convert string "NaN" to MATLAB NaN
if isfield(settings, 'Threshold') && strcmp(settings.Threshold, 'NaN')
    settings.Threshold = NaN;
end

% Create a full identifier string based on the options
% fullIdentifier = createIdentifier(settings);
end