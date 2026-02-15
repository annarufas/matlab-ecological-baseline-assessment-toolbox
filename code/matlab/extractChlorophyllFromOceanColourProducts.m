function [S,prodLabel,prodName,nChlaProducts] = ...
    extractChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci)

% EXTRACTCHLOROPHYLLFROMOCEANCOLOURPRODUCTS Extract chlorophyll data from 
% three ocean colour products (CMEMS, NASA, and OC-CCI) and organise the
% data into a structure for easy access and analysis.
%
%   INPUT: 
%       AScmems       - table with time-series data from CMEMS for our area of study
%       ASnasa        - table with time-series data from NASA for our area of study
%       ASoccci       - table with time-series data from OC-CCI for our area of study
%
%   OUTPUT:
%       S             - data structure with chlorophyll data extracted from selected OC datasets
%       prodLabel     - string containing selected OC product labels for plotting
%       prodName      - string containing selected OC product names
%       nChlaProducts - number of products extracted
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 4 May 2024
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Presets

% Define product labels and names for selected ocean colour products
prodLabel = {'4 km, global, OLCI',...
             '4 km, global, Aqua-MODIS',...
             '4 km, global, VIIRS-SNPP',...
             '4 km, global, OC-CCI multi',...
             '1 km, global, OC-CCI multi',...
             '1 km, regional, CMEMS multi',...
             '300 m, regional, OLCI',...
             '7 km, regional, reanalysis'};

prodName = {'global_olci_r4km',...
            'global_aquamodis_r4km',...
            'global_viirssnpp_r4km',...
            'global_occcimulti_r4km',...
            'global_occcimulti_r1km',...
            'regional_cmemsmulti_r1km',...
            'regional_olci_r300m',...
            'regional_cmemsreanalysis_r7km'};

nChlaProducts = 8;

% Initialise a structure to hold the extracted chlorophyll data
S = struct();

%% Extract chlorophyll data from products

% Retrieve ID field names and dataset names
[cmemsLoadedDatasetNames] = getTimeSeriesDatasetNames(AScmems);
[nasaLoadedDatasetNames] = getTimeSeriesDatasetNames(ASnasa);
[occciLoadedDatasetNames] = getTimeSeriesDatasetNames(ASoccci);
     
for iDataset = 1:numel(cmemsLoadedDatasetNames)
    switch cmemsLoadedDatasetNames{iDataset}
        case 'obs_satell_glob_cmems_olci_4km_plk'
            idxProd = find(strcmp({AScmems.ID}, 'obs_satell_glob_cmems_olci_4km_plk'));
            S.global_olci_r4km.dataset = AScmems(idxProd).dataset(:,:,:,1);
            S.global_olci_r4km.time = AScmems(idxProd).time; 
            S.global_olci_r4km.lat = AScmems(idxProd).lat; 
            S.global_olci_r4km.lon = AScmems(idxProd).lon; 
        case 'obs_satell_reg_cmems_multi_1km_plk'
            idxProd = find(strcmp({AScmems.ID}, 'obs_satell_reg_cmems_multi_1km_plk'));
            S.regional_cmemsmulti_r1km.dataset = AScmems(idxProd).dataset(:,:,:,1);
            S.regional_cmemsmulti_r1km.time = AScmems(idxProd).time;
            S.regional_cmemsmulti_r1km.lat = AScmems(idxProd).lat;
            S.regional_cmemsmulti_r1km.lon = AScmems(idxProd).lon;
        case 'obs_satell_reg_cmems_olci_300m_plk'
            idxProd = find(strcmp({AScmems.ID}, 'obs_satell_reg_cmems_olci_300m_plk'));
            S.regional_olci_r300m.dataset = AScmems(idxProd).dataset(:,:,:,1);
            S.regional_olci_r300m.time = AScmems(idxProd).time; 
            S.regional_olci_r300m.lat = AScmems(idxProd).lat;
            S.regional_olci_r300m.lon = AScmems(idxProd).lon;
        case 'mod_bgc_reg_chl'
            idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_chl'));
            S.regional_cmemsreanalysis_r7km.dataset = AScmems(idxProd).dataset(:,:,:,1);
            S.regional_cmemsreanalysis_r7km.time = AScmems(idxProd).time; 
            S.regional_cmemsreanalysis_r7km.lat = AScmems(idxProd).lat;
            S.regional_cmemsreanalysis_r7km.lon = AScmems(idxProd).lon;
    end
end

for iDataset = 1:numel(nasaLoadedDatasetNames)
    switch nasaLoadedDatasetNames{iDataset}
        case 'aquamodis_4km'
            idxProd = find(strcmp({ASnasa.ID}, 'aquamodis_4km'));
            S.global_aquamodis_r4km.dataset = ASnasa(idxProd).dataset(:,:,:);
            S.global_aquamodis_r4km.time = ASnasa(idxProd).time;
            S.global_aquamodis_r4km.lat = ASnasa(idxProd).lat; 
            S.global_aquamodis_r4km.lon = ASnasa(idxProd).lon; 
        case 'viirssnpp_4km'
            idxProd = find(strcmp({ASnasa.ID}, 'viirssnpp_4km'));
            S.global_viirssnpp_r4km.dataset = ASnasa(idxProd).dataset(:,:,:);
            S.global_viirssnpp_r4km.time = ASnasa(idxProd).time;
            S.global_viirssnpp_r4km.lat = ASnasa(idxProd).lat;
            S.global_viirssnpp_r4km.lon = ASnasa(idxProd).lon;
    end
end

for iDataset = 1:numel(occciLoadedDatasetNames)
    switch occciLoadedDatasetNames{iDataset}
        case 'occci_1km_1day'
            idxProd = find(strcmp({ASoccci.ID}, 'occci_1km_1day'));
            S.global_occcimulti_r1km.dataset = ASoccci(idxProd).dataset(:,:,:);
            S.global_occcimulti_r1km.time = ASoccci(idxProd).time; 
            S.global_occcimulti_r1km.lat = ASoccci(idxProd).lat; 
            S.global_occcimulti_r1km.lon = ASoccci(idxProd).lon; 
        case 'occci_4km_1day'
            idxProd = find(strcmp({ASoccci.ID}, 'occci_4km_1day'));
            S.global_occcimulti_r4km.dataset = ASoccci(idxProd).dataset(:,:,:);
            S.global_occcimulti_r4km.time = ASoccci(idxProd).time; 
            S.global_occcimulti_r4km.lat = ASoccci(idxProd).lat; 
            S.global_occcimulti_r4km.lon = ASoccci(idxProd).lon; 
    end
end

end % extractChlorophyllFromOceanColourProducts