function [HPLC] = prepareCefasHplcData(filenameCefasHplcRawData,...
    AVG_EUPHOTIC_ZONE_DEPTH,AVG_SHELF_DEPTH,pathBathymetryFile)

% PREPARECEFASHPLCDATA Create a working file from the CEFAS' raw HPLC dataset.
% For that, remove observations at depths deeper than AVG_EUPHOTIC_ZONE_DEPTH 
% as well as those in a water column with a bathymetry deeper AVG_SHELF_DEPTH. 
% If there are multiple observations within the first AVG_EUPHOTIC_ZONE_DEPTH m 
% of the water column for a location and time, only retain the observation 
% closest to the surface.
%
%   INPUT: 
%       filenameCefasHplcRawData - .csv file containing the raw HPLC data
%       AVG_EUPHOTIC_ZONE_DEPTH  - depth in m
%       AVG_SHELF_DEPTH          - depth in m
%       pathBathymetryFile       - .nc file containing the ETOPO bathymetric dataset
%
%   OUTPUT:
%       HPLC                     - Matlab table with filtered HPLC data
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 24 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

% Abbreviations used in the dataset:
% ug/L              = mg/m3
% TP                total pigments
% TChl              total chlorophyll (TChla + Chlb + Chlc)
% TChla             total chlorophyll a
% AP                accessory pigments (TC + Chlb + Chlc)
% TC                total carotenoids (Allo + But + Cara + Carb + Carab + Diad + Diato + Fuc + Hex + Lut + etc.)
% Degrad_products   Pheidea + Pheoa

fprintf('\nFiltering the HPLC dataset...\n')

%% Read the raw CEFAS HPLC dataset

% Import the datetime stamp as a string and then apply some formatting. 
% Notice that datetime comes in a mix of formats (sometimes with AM/PM, 
% sometimes not)
T = readtable(fullfile('.','data','raw','CEFAS_HPLC',filenameCefasHplcRawData),'DatetimeType','text');
T.DateTime = datetime(T.DateTime,'InputFormat','MM/dd/yyyy h:mm:ss a');
cefasHplc = T;
% Add ID column
T2 = addvars(T,(1:height(T))','Before','Survey_name','NewVariableNames','idd');
cefasHplc = T2;

%% Process it

idxColStartPigmentVars = 9;

% Find unique combinations of latitude x longitude x depth x time. 
% The ids that are non-unique (i.e., repeated) are replicates, most of them
% flagged with "NFS" (CEFAS does not specify meaning of that flag).
[C,ia,ic] = unique([cefasHplc.Latitude cefasHplc.Longitude... 
    cefasHplc.Sample_depth string(cefasHplc.DateTime)],'rows');

% Compute the average of the replicates
avgReplicate = table();
iColNew = 0;
for iColOld = idxColStartPigmentVars:width(cefasHplc)
    iColNew = iColNew + 1;
    avgReplicate.(iColNew) = accumarray(ic, cefasHplc.(iColOld), [], @mean);
end

% Take variable names for pigment variables
avgReplicate.Properties.VariableNames = cefasHplc.Properties.VariableNames(idxColStartPigmentVars:end); 

% Extract structural variables
structVars = cefasHplc(ia,1:(idxColStartPigmentVars-1)); % idxColStartPigmentVars-1 = 8

% Combine structural + numerical sections of the table
hplcNoRep = [structVars avgReplicate];

%% Only take surface values

% Only take those locations where sampling was conducted in the first
% surface depth metres
idxsSurfaceDepth = find(hplcNoRep.Sample_depth <= AVG_EUPHOTIC_ZONE_DEPTH);
hplcNoRepSurf = hplcNoRep(idxsSurfaceDepth,:);

% When there's more than one observation in the first surface depth metres 
% for a specific location and time, take the observation closest to the surface
[C,ia,ic] = unique([hplcNoRepSurf.Latitude hplcNoRepSurf.Longitude... 
    string(hplcNoRepSurf.DateTime)],'first','rows');
hplcNoRepSurfUnique = hplcNoRepSurf(ia,:);

%% Remove locations that are sitting in deep water columns or on land

[HPLC] = addBathymetryAndFilter(hplcNoRepSurfUnique,'Longitude','Latitude',...
    pathBathymetryFile,AVG_SHELF_DEPTH);

%% Add season 

HPLC = addSeason(HPLC, 'DateTime');

%% Save the data

save(fullfile('.','data','processed','cefasHPLCfiltered.mat'),'HPLC')

% Save data to a .csv file that can be read into Python to extract data for matchups
%csvTable = P(:,[1,5,6,7,9,12]); % 1=idd, 5=DateTime, 6=Latitude, 7=Longitude, 9=Sample_depth, 12=TChlA_ug_L
csvTable = HPLC;
writetable(csvTable,fullfile('.','data','processed','cefasHPLCfiltered.csv'),'Delimiter','comma')

fprintf('\n... done.\n')

end % prepareCefasHplcData
