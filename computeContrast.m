%Calculates the sum envelope contrast of range profiles
function Contrast_C = computeContrast(HRRProfiles)
    SumEnvelope_S = sum(abs(HRRProfiles),1);
    Contrast_C = sum(SumEnvelope_S.*SumEnvelope_S);
end