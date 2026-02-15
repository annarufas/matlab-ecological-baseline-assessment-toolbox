function [rmse,log_rmse] = calcRootMeanSquaredError(S, I)

N = length(S);

squaredErrors = (S-I).^2;
sumOfSquaredErrors = sum(squaredErrors);
rmse = sqrt(sumOfSquaredErrors/N);

squaredLogErrors = (log(S)-log(I)).^2; 
sumOfSquaredLogErrors = sum(squaredLogErrors);
log_rmse = sqrt(sumOfSquaredLogErrors/N);

end