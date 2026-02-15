function [regressionStatistics] = calculateMatchupRegressionStatistics(...
    insituData,estimatedData)

% CALCULATEMATCHUPREGRESSIONSTATISTICS Puts together a series of regression 
% statistics for matchup analysis.
%
%   INPUT:
%       insituData           - in situ data
%       estimatedData        - estimated data
%
%   OUTPUT:
%       regressionStatistics - regressionStatistics
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 25 April 2024   
%
% =========================================================================
%%

N = length(estimatedData);

% Fit line
regressCoeff = polyfit(insituData,estimatedData,1);
slope = regressCoeff(1);
intercept = regressCoeff(2);

% Correlation explains the strength of the relationship between an independent and dependent variable
[corrCoeff,pvalue] = corr(insituData,estimatedData,'Type','Pearson','Rows','complete'); % Pearson linear correlation coefficient using only the rows that contain no missing values.

% R-squared explains to what extent the variance of one variable explains the variance of the second variable
rsquared = corrCoeff^2;

[rmse,log_rmse] = calcRootMeanSquaredError(estimatedData,insituData); % greek letter: phi
[me,log_me] = calcMeanError(estimatedData,insituData); % aka bias, greek letter: delta
[mae,log_mae] = calcMeanAbsoluteError(estimatedData,insituData);
mrpe = calcMeanRelativePercentageError(estimatedData,insituData);
mape = calcMeanAbsolutePercentageError(estimatedData,insituData);

% Add regresion statistics to the table
regressionStatistics = [N,slope,intercept,rsquared,corrCoeff,rmse,log_rmse,me,log_me,mae,log_mae,mrpe,mape];

end