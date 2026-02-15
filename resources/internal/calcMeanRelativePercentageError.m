function mrpe = calcMeanRelativePercentageError(S, I)

N = length(I);
relativePercentageError = 100.*((S-I)./I);
sumOfRelativePercentageError = sum(relativePercentageError);
mrpe = sumOfRelativePercentageError/N;

end
