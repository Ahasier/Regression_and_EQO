function useExtraFeatures = useExtraFeatures(settings)
if isfield(settings, 'usePhylogeny') && strcmp(settings.usePhylogeny, 'On')
    useExtraFeatures = true;
else
    useExtraFeatures = false;
end
end