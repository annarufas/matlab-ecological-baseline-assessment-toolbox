function createAreaStudyShapefile(radiusAreaStudy,...
    fullPathUkOffshoreGcsLicensesShapefile,fullPathGcsSiteShapefileDir,...
    fullPathAreaStudyShapefileDir,fileNameGcsSiteCentralCoords,...
    shapefileNameGcsSite,shapefileNameAreaStudy)

% CREATEAREASTUDYSHAPEFILE Create a shapefile (.shp) for a GCS site and its 
% monitoring footprint (area of study).
%
%   INPUT: 
%       radiusAreaStudy                        - radius of the area of study in kilometers
%       fullPathUkOffshoreGcsLicensesShapefile - path to shapefile containing NSTA's UKCS offshore licensed sites
%       fullPathGcsSiteShapefileDir            - directory to store our GCS site shapefile
%       fullPathAreaStudyShapefileDir          - directory to store our area of study shapefile
%       fileNameGcsSiteCentralCoords           - .mat file with our GCS site central coordinates
%       shapefileNameGcsSite                   - .shp file with our GCS site coordinates
%       shapefileNameAreaStudy                 - .shp file with our area of study coordinates
%
%   OUTPUT:
%       Creates two shapefiles:
%       1. A shapefile for the specified GCS site.
%       2. A shapefile for the area of study around the GCS site, based on the radius.
%
%   This script uses this external function:
%       m_shaperead - from m_map, FileExchange

%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 19 April 2024
%   Version 1.1 - Updated 6 Jan 2025: simplified terminology
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Create a shapefile for the GCS site Endurance and calculate central coordinates

if ~exist(fullfile(fullPathGcsSiteShapefileDir,shapefileNameGcsSite),'file')
    
    fprintf('\nCreating shapefile for our GCS site...')
    
    % Change accordingly if not using GCS site Endurance (Endurance is the first site in the shapefile in fullPathUkOffshoreGcsLicensesCoords)
    idxEndurance = 1; 

    % Read UK Continental Shelf Carbon Storage sites
    ukcsGcsSites = shaperead(fullPathUkOffshoreGcsLicensesShapefile); % lat/lon coordinates

    % Extract coordinates for our GCS site
    gcsSiteLon = ukcsGcsSites.ncst{idxEndurance}(:,1);  
    gcsSiteLat = ukcsGcsSites.ncst{idxEndurance}(:,2);
    
    % Save the GCS site shapefile
    gcsSite.Geometry = 'Polygon';
    gcsSite.X = gcsSiteLon; 
    gcsSite.Y = gcsSiteLat;  
    gcsSite.Name = 'Rectangle'; 
    shapewrite(gcsSite,fullfile(fullPathGcsSiteShapefileDir,shapefileNameGcsSite,'.shp'))

    % Calculate central coordinates for our GCS site
    pgon = polyshape(gcsSiteLon,gcsSiteLat);
    [gcsSiteLonCentre,gcsSiteLatCentre] = centroid(pgon);
    
    % Save central coordinates as .mat file
    save(fullfile(fullPathGcsSiteShapefileDir,fileNameGcsSiteCentralCoords),...
        'gcsSiteLonCentre','gcsSiteLatCentre')
    
    fprintf('...done.\n')
    
end

%% Create a shapefile for the broader area of study around the GCS site Endurance
   
% Create a boundary box region expanding radiusAreaStudy km to the north 
% and south and radiusAreaStudy km to the east and west of the GCS site. 
% This function creates a shapefile with the lat/lon coordinates 
% corrresponding to the RADIUS_AREA_STUDY km expansion centred around the 
% GCS site.

if ~exist(fullfile(fullPathAreaStudyShapefileDir,shapefileNameAreaStudy),'file')
    
    fprintf('Creating shapefile for our broader area of study...')

    load(fullfile(fullPathGcsSiteShapefileDir,fileNameGcsSiteCentralCoords),...
        'gcsSiteLonCentre','gcsSiteLatCentre')
   
    ydistDeg_initguess = 2; % how many km are 2ºN?
    xdistDeg_initguess = 2; % how many km are 2ºE?
    distStepDeg = 0.001;    % step distance in degrees

    % The following uses an iteration method to approximate the distance to 
    % radiusAreStudy

    % First, minimise xdistDeg
    [arclen,~] = distance([gcsSiteLatCentre, gcsSiteLonCentre],... 
        [gcsSiteLatCentre, gcsSiteLonCentre+xdistDeg_initguess]);
    xdistDeg = xdistDeg_initguess;
    distkm = deg2km(arclen);
    while distkm > radiusAreaStudy
        xdistDeg = xdistDeg - distStepDeg;
        [arclen,~] = distance([gcsSiteLatCentre, gcsSiteLonCentre],... 
            [gcsSiteLatCentre, gcsSiteLonCentre+xdistDeg]);
        distkm = deg2km(arclen);
        if distkm <= radiusAreaStudy
            break
        end
    end

    fprintf('\nThe best minimised distance in the x direction is %5.3f km, which corresponds to %5.3fºE.',distkm,xdistDeg)

    % Second, minimise ydistDeg
    [arclen,~] = distance([gcsSiteLatCentre, gcsSiteLonCentre],... 
        [gcsSiteLatCentre+ydistDeg_initguess, gcsSiteLonCentre]);
    ydistDeg = ydistDeg_initguess;
    distkm = deg2km(arclen);
    while distkm > radiusAreaStudy
        ydistDeg = ydistDeg - distStepDeg;
        [arclen,~] = distance([gcsSiteLatCentre, gcsSiteLonCentre],... 
            [gcsSiteLatCentre+ydistDeg, gcsSiteLonCentre]);
        distkm = deg2km(arclen);
        if distkm <= radiusAreaStudy
            break
        end
    end

    fprintf('\nThe best minimised distance in the y direction is %5.3f km, which corresponds to %5.3fºN.',distkm,ydistDeg)

    bboxlat = zeros(5,1);
    bboxlon = zeros(5,1);

    bboxlat(1) = gcsSiteLatCentre - ydistDeg;
    bboxlon(1) = gcsSiteLonCentre - xdistDeg;
    bboxlat(2) = gcsSiteLatCentre + ydistDeg;
    bboxlon(2) = bboxlon(1);
    bboxlat(3) = bboxlat(2);
    bboxlon(3) = gcsSiteLonCentre + xdistDeg;
    bboxlat(4) = bboxlat(1);
    bboxlon(4) = bboxlon(3);
    bboxlat(5) = bboxlat(1);
    bboxlon(5) = bboxlon(1);

    bbox.Geometry = 'Polygon';
    bbox.X = bboxlon; 
    bbox.Y = bboxlat;  
    bbox.Name = 'Rectangle'; 
    shapewrite(bbox,fullfile(fullPathAreaStudyShapefileDir,strcat(shapefileNameAreaStudy,'.shp')))

    fprintf('\n...done.\n')

end

end % createAreaStudyShapefile
