function [sceneAverageNoOutliers,sceneNumPixels] = calcSceneAverage(dataset,time)

% PLOTHOVMOLLERDIAGRAMSFORDRIVERS Calculate scene average, remove oulier
% data and calculate monthly and yearly averages.
%
%   INPUT: 
%       dataset           - table with time-series data from CMEMS for our area of study
%       time
%
%   OUTPUT:
%       sceneAverageNoOutliers  - data structure with data extracted from selected datasets
%       sceneNumPixels          - string containing selected product labels for plotting
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 7 May 2024
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Scene average

sceneAverage = NaN(numel(time),1);
sceneNumPixels = NaN(numel(time),1);
for iTimeStep = 1:numel(time)
    sceneAverage(iTimeStep) = mean(dataset(:,:,iTimeStep),'all','omitnan');
    sceneNumPixels(iTimeStep) = nnz(~isnan(dataset(:,:,iTimeStep)));            
end

%% Consider removal of outliers using a moving median of 7 days

sceneAverageNoOutliers = filloutliers(sceneAverage,'clip','movmedian',days(7),'SamplePoints',time);

end % calcSceneAverage