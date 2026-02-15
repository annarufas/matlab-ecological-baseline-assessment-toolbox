function [PH] = prepareInsituPhData(filenameSsbProgramPhData,filenameUkoaProgramPhData,...
    AVG_SHELF_DEPTH,pathBathymetryFile)

% PREPAREINSITUPHDATA Read in ALK and DIC data from the SSB and UKOA
% programs, calculate pH using Matlab's CO2SYS tool and filter out data
% points that are not sitting on the shelf seas.
%
%   INPUT:
%       filenameSsbProgramPhData  - .xlsx file with ALK and DIC data from the SSB program
%       filenameUkoaProgramPhData - .xlsx file with ALK and DIC data from the UKOA program
%       AVG_SHELF_DEPTH           - depth in m
%       pathBathymetryFile        - .nc file containing the ETOPO bathymetric dataset
%
%   OUTPUT:
%       PH                        - Matlab table with combined data from both programs and filtered pH data
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

fprintf('\nReading in ALK and DIC data from the SSB and UKOA programs,')
fprintf('\ncalculating pH using CO2SYS and filtering the data...\n')

%% Read in the .xlsx files into tables and do some formatting

T_ssb = readtable(fullfile('.','data','raw','Greenwood_pH',filenameSsbProgramPhData));
T_ukoa = readtable(fullfile('.','data','raw','Greenwood_pH',filenameUkoaProgramPhData));

% Add tag column
T_ssb = addvars(T_ssb,repmat({'ssb'},[height(T_ssb) 1]),'Before','Agency','NewVariableNames','Program');
T_ukoa = addvars(T_ukoa,repmat({'ukoa'},[height(T_ukoa) 1]),'Before','ID','NewVariableNames','Program');

% For T_ukoa, calculate the mean of TA and DIC from the three replicates
taColumnNames = {'TA_rep1_umol_kg','TA_rep2_umol_kg','TA_rep3_umol_kg'};
taMean = mean(T_ukoa{:,taColumnNames},2,'omitnan');
dicColumnNames = {'DIC_rep1_umol_kg','DIC_rep2_umol_kg','DIC_rep3_umol_kg'};
dicMean = mean(T_ukoa{:,dicColumnNames},2,'omitnan');
T_ukoa = addvars(T_ukoa,taMean,'Before','DIC_rep1_umol_kg','NewVariableNames','TA_mean_umol_kg');
T_ukoa = addvars(T_ukoa,dicMean,'Before','SPM_mg_L','NewVariableNames','DIC_mean_umol_kg'); 

%% Calculate pH using CO2SYS from DIC and TAlk

talk        = [T_ssb.TA_umol_kg;T_ukoa.TA_mean_umol_kg]; % umol/kg
tdic        = [T_ssb.DIC_umol_kg;T_ukoa.DIC_mean_umol_kg]; % umol/kg
salinity    = [T_ssb.Salinity_PSU;T_ukoa.Salinity_PSU]; % PSU
temperature = [T_ssb.Temperature_degC;T_ukoa.Temperature_degC]; % degC
pressure    = [T_ssb.Sample_depth_m;T_ukoa.Pressure_dbar]; % dbar
silicate    = [T_ssb.Silicate_umol_kg;T_ukoa.Si_umol_kg]; % umol/kg
phosphate   = [T_ssb.Phosphate_umol_kg;T_ukoa.PO4_umol_kg]; % umol/kg

par1type = 1; % The first parameter supplied is of type "1", which is "alkalinity"
par2type = 2; % The first parameter supplied is of type "2", which is "DIC"
pHscale  = 1; % pH scale at which the input pH is reported ("1" means "Total Scale" and it is the default)
k1k2c    = 4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = 1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

pH = NaN(size(talk));
for i = 1:height(talk)   
    [A,headers] = CO2SYS(talk(i),tdic(i),par1type,par2type,salinity(i),...
        temperature(i),temperature(i),pressure(i),pressure(i),silicate(i),...
        phosphate(i),pHscale,k1k2c,kso4c);
    pH(i) = A(18);
end

% Add the pH into the tables
T_ssb = addvars(T_ssb,pH(1:height(T_ssb)),'Before','Silicate_umol_kg','NewVariableNames','pH');
T_ukoa = addvars(T_ukoa,pH(height(T_ssb)+1:end),'Before','SPM_mg_L','NewVariableNames','pH'); 

%% Average replicates

ssbPhColIdx = find(strcmp(T_ssb.Properties.VariableNames,'pH'));
ukoaPhColIdx = find(strcmp(T_ukoa.Properties.VariableNames,'pH'));
ssbIdxStructVars = [1,3,4,7,8,9];
ukoaIdxStructVars = [1,6,8,9,10,13];

% Process T_ssb dataset
[ssbPhAvgReplicates, ssbStructVars] = aggregateReplicates(T_ssb,...
    ssbPhColIdx, ssbIdxStructVars);

% Combine structural + numerical sections of the table
ssbPhNoRep = [ssbStructVars, ssbPhAvgReplicates];

% Process T_ukoa dataset
[ukoaPhAvgReplicates, ukoaStructVars] = aggregateReplicates(T_ukoa,...
    ukoaPhColIdx, ukoaIdxStructVars);

% Combine structural + numerical sections of the table
ukoaPhNoRep = [ukoaStructVars, ukoaPhAvgReplicates];

% Combine tables
insituPhNorthSea = [ssbPhNoRep; ukoaPhNoRep];

% Delet rows where pH is NaN
isPhNaN = isnan(insituPhNorthSea.pH);
rowsToKeep = ~isPhNaN;
insituPhNorthSea = insituPhNorthSea(rowsToKeep, :);

%% Remove locations that are sitting in deep water columns or on land

[PH] = addBathymetryAndFilter(insituPhNorthSea,'Longitude_degE','Latitude_degN',...
    pathBathymetryFile,AVG_SHELF_DEPTH);

%% Add season 

PH = addSeason(PH,'DateTime');

%% Save the data

save(fullfile('.','data','processed','greenwoodPhFiltered.mat'),'PH')

% Save data to a .csv file that can be read into Python to extract data for matchups
csvTable = PH;
writetable(csvTable,fullfile('.','data','processed','greenwoodPhFiltered.csv'),'Delimiter','comma')

fprintf('\n... done.\n')

%%
% -------------------------------------------------------------------------
% LOCAL FUNCTION
% -------------------------------------------------------------------------

% *************************************************************************

function [avgReplicates, structVars] = aggregateReplicates(dataTable,...
    phColIdx,structVarIndices)

    % Find unique combinations of latitude x longitude x depth x time.
    [C,ia,ic] = unique([dataTable.Latitude_degN dataTable.Longitude_degE...
        dataTable.Sample_depth_m string(dataTable.DateTime)], 'rows');
    
    % Compute the average of the replicates
    avgReplicates = accumarray(ic, dataTable.(phColIdx), [], @mean);
    avgReplicates = table(avgReplicates, 'VariableNames', {'pH'});
    
    % Extract structural variables
    structVars = dataTable(ia, structVarIndices);
    
end

% *************************************************************************

end % prepareInsituPhData
