function [bicepData] = addDepthIntegrationProductsToBicepArray(bicepData,...
  pathMldDir,pathAreaStudyShapefile,cmemsData,filenameBicepDataProcessed)

% ADDDEPTHINTEGRATIONPRODUCTSTOBICEPARRAY Integrate the mixed layer depth
% (MLD) product from MIMOC and the euphotic zone depth (Zeu) product from
% CMEMS into the BICEP data structure. These products are necessary for
% calculating carbon-integrated quantities from the BICEP's products.
%
%   INPUT:
%       bicepData                  - Matlab structure with BICEP data
%       pathMldDir                 - path to MLD file
%       cmemsData                  - Matlab structure with CMEMS data
%       pathAreaStudyShapefile     - shapefile with our area of study
%       filenameBicepDataProcessed - .mat file containing bicepData
%
%   OUTPUT:
%       bicepData  - Matlab table with BICEP data, now containing MIMOC's MLD and CMEMS's Zeu
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 3 May 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Extract MLD data from MIMOC (monthly climatology)
 
fprintf("\nReading MLD from MIMOC and adding it to the BICEP data array...")

fileNames = dir(fullfile(pathMldDir,'*.nc'));

% Initialise latitude and longitude only once outside the month loop
monthFolderPath = fullfile(pathMldDir, fileNames(1).name);
lat = ncread(monthFolderPath, 'LATITUDE'); % -89 to 89ºN
lon = ncread(monthFolderPath, 'LONGITUDE'); % 0-360º

% Transform longitude to degrees East
lon = mod(lon + 180, 360) - 180; % -180 to 180ºE

% Adjust coordinates for the area of study
[idxMinLat,idxMaxLat,idxMinLon,idxMaxLon,lat,lon] =...
    adjustAreaStudyCoordinates(lat,lon,pathAreaStudyShapefile);

% Initialise the output data array to store 12 months of data
Dout = NaN(numel(lat),numel(lon),12);
        
for iMonth = 1:12
    
    monthFolderPath = fullfile(pathMldDir, fileNames(iMonth).name);
    %S = ncinfo(monthFolderPath); 
    
    % Read data array
    Dtemp = ncread(monthFolderPath,'DEPTH_MIXED_LAYER',...
        [idxMinLon, idxMinLat], [idxMaxLon-idxMinLon+1, idxMaxLat-idxMinLat+1]); 
    
    % Permute dimensions (swap lat and lon) 
    Dperm = permute(Dtemp, [2, 1]); 

    % Store the sorted data in the output array
    Dout(:,:,iMonth) = Dperm;
    
end

% Find latest position occupied in the structure carbonStocks
populatedIDs = ~cellfun(@isempty, {bicepData.ID});
iLatestPositionOccupied = sum(populatedIDs);
iPositionFree = iLatestPositionOccupied + 1;

% Save information into output array
bicepData(iPositionFree).ID = 'mld_mimoc_poc';
bicepData(iPositionFree).varNames = 'mld';
bicepData(iPositionFree).units = 'm';
bicepData(iPositionFree).dataset = Dout;
bicepData(iPositionFree).lat = double(lat);
bicepData(iPositionFree).lon = double(lon);

fprintf("\n ... done.\n")

%% Extract Zeu from CMEMS

disp("\nReading Zeu from CMEMS and adding it to the BICEP data array...")

[cmemsZeu] = calculateZeuFromCmemsKd(cmemsData);

% Find latest position occupied in the structure carbonStocks
populatedIDs = ~cellfun(@isempty, {bicepData.ID});
iLatestPositionOccupied = sum(populatedIDs);
iPositionFree = iLatestPositionOccupied + 1;

% Save information into output array
bicepData(iPositionFree).ID = 'cmems_zeu';
bicepData(iPositionFree).varNames = 'zeu';
bicepData(iPositionFree).units = 'm';
bicepData(iPositionFree).dataset = cmemsZeu;
bicepData(iPositionFree).lat = double(cmemsData(idxKdCmems).lat);
bicepData(iPositionFree).lon = double(cmemsData(idxKdCmems).lon);
bicepData(iPositionFree).time = cmemsData(idxKdCmems).time;

fprintf("\n... done.")

%% Save outputs

fprintf("\n... finished creating BICEP data structure, saving it...")
save(fullfile('.','data','processed',filenameBicepDataProcessed),'bicepData','-v7.3')
fprintf("\n... saving completed.\n")

end % addDepthIntegrationProductsToBicepArray