function isReal = usingRealData(settings)
isReal = isfield(settings, 'RealData') && strcmp(settings.RealData, 'On');
end