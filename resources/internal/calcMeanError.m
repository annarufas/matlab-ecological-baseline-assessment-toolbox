function [me,log_me] = calcMeanError(S, I)

% Also called mean difference (MD), mean error (ME), mean bias (MB) or just
% bias
% Greek symbol: delta

N = length(I);

me = sum(S-I)/N; 

log_me = sum(log(S)-log(I))/N;

end