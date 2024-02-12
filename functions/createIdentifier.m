function fullIdentifier = createIdentifier(settings)
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