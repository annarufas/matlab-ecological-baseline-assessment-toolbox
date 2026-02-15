function [mae,log_mae] = calcMeanAbsoluteError(S, I)

% Also called mean absolute deviation or difference (MAD) or error (MAE)

N = length(I);

mae = sum(abs(S-I))/N; 

log_mae = sum(abs(log(S)-log(I)))/N;

end