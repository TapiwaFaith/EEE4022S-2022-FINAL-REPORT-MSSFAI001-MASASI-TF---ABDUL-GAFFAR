%Calculates the sum envelope contrast of range profiles
function Entropy_H = computeEntropy(HRRProfiles)
    SumEnvelope_S = sum(abs(HRRProfiles),1);
    Entropy_H = -(sum(SumEnvelope_S.*log(SumEnvelope_S)));
end