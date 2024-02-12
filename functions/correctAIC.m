function aicc = correctAIC(aic, numParameters, numSamples)
aicc = aic + 2*(numParameters*(numParameters + 1))/(numSamples - numParameters - 1);
end