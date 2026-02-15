function [idxMinLat,idxMaxLat,idxMinLon,idxMaxLon,latVector,lonVector] =... 
    adjustAreaStudyCoordinates(latVector,lonVector,pathAreaStudyShapefile)

% ADJUSTAREASTUDYCOORDINATES Takes in latitude and longitude vectors as 
% and calculates the indices of the closest values to the bounding box 
% defined by the shapefile.
%
%   INPUT:
%       latVector              - latitude vector (degress north)
%       lonVector              - longitude vector (degrees east)
%       pathAreaStudyShapefile - shapefile with our area of study
% 
%   OUTPUT:
%       idxMinLat              - index corresponding to min latitude in area of study
%       idxMaxLat              - index corresponding to max latitude in area of study
%       idxMinLon              - index corresponding to min longitude in area of study
%       idxMaxLon              - index corresponding to max longitude in area of study
%       latVector              - adjusted latitude vector (degress north)
%       lonVector              - adjusted longitude vector (degrees east)
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 28 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Box of 100 x 100 km around the Endurance site

areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates
minLatAreaStudy = areaStudy.MBRy(1); % 53.758
maxLatAreaStudy = areaStudy.MBRy(2); % 54.656
minLonAreaStudy = areaStudy.MBRx(1); % 0.265
maxLonAreaStudy = areaStudy.MBRx(2); % 1.801

%% Find closest index in the coordinate vectors

[~, closestIdxMinLat] = min(abs(latVector - minLatAreaStudy));
[~, closestIdxMaxLat] = min(abs(latVector - maxLatAreaStudy));
[~, closestIdxMinLon] = min(abs(lonVector - minLonAreaStudy));
[~, closestIdxMaxLon] = min(abs(lonVector - maxLonAreaStudy));

%% Expand the ranges around the closest indices 
% This ensures a buffer around the area of study. Additionally, account for 
% the fact that latitude is often unsorted (longitude is not)

if (issorted(latVector))
    idxMinLat = closestIdxMinLat-1;
    idxMaxLat = closestIdxMaxLat+1;  
    latVector = latVector(idxMinLat:idxMaxLat);
elseif (~issorted(latVector))
    idxMinLat = closestIdxMinLat+1;
    idxMaxLat = closestIdxMaxLat-1;
    latVector = latVector(idxMaxLat:idxMinLat);
end
idxMinLon = closestIdxMinLon-1;
idxMaxLon = closestIdxMaxLon+1;
lonVector = lonVector(idxMinLon:idxMaxLon);
  
end % adjustAreaStudyCoordinates 
