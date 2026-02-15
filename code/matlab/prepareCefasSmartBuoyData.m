function [SB] = prepareCefasSmartBuoyData(filenamesCefasSmartBuoyData)

% PREPARECEFASSMARTBUOYDATA Reads in the SmartBuoy dataset.
%
%   INPUT: 
%       filenamesCefasSmartBuoyData - string containing the names of three files with SmartBuoy data
%
%   OUTPUT:
%       SB – table with SmartBuoy data
%         
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 19 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Import data into a table

nFiles = numel(filenamesCefasSmartBuoyData);

varsToImport = {'par',...
                'deployment_group',...
                'dateTime',...
                'lat',...
                'lon',...
                'depth',...
                'value',...
                'stdev',...
                'n',...
                'unit'};

renameVarsToImport = {'variable',...
                      'buoy',...
                      'dateTime',...
                      'lat',...
                      'lon',...
                      'depth',...
                      'value',...
                      'std',...
                      'N',...
                      'units'};

SB = table();

for iFile = 1:nFiles

    thisFile = fullfile('.','data','raw','CEFAS_SmartBuoy',filenamesCefasSmartBuoyData{iFile});

    % Only import specific columns and make sure that they will be
    % imported with the same variable type
    opts = detectImportOptions(thisFile);
    opts.SelectedVariableNames = varsToImport;
    opts = setvartype(opts,varsToImport,...
        {'categorical','categorical','datetime','double','double','double',...
        'double','double','double','categorical'});

    SBnew = readtable(thisFile,opts); 
    SBnew = renamevars(SBnew,varsToImport,renameVarsToImport);
    SB = [SB; SBnew];

end

%% Add season 

SB = addSeason(HPLC, 'dateTime');

%% Save the dataset

save(fullfile('.','data','processed','cefasSmartBuoy.mat'),'SB')

end % prepareCefasSmartBuoyData
