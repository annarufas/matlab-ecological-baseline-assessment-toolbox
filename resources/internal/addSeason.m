function dataTable = addSeason(dataTable, dateTimeColName)

% ADDSEASON Assign season to datetime values.
%
%   INPUT:
%       dataTable       - table containing datetime values
%       dateTimeColName - name of the column with datetime values
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

%% Initialise the 'season' column as a categorical array of NaN

dataTable.season = categorical(repelem(NaN, height(dataTable), 1));

%% Definitions

% Define the start and end years based on the datetime column
years = year(dataTable.(dateTimeColName));
startYear = min(years);
endYear = max(years);

% Define the start and end dates for each season
seasons = {'Winter', 'Spring', 'Summer', 'Autumn'};
seasonStartEnd = {
    datetime(startYear-1,12,21), datetime(startYear,3,20); % Winter
    datetime(startYear,3,21),    datetime(startYear,6,20); % Spring
    datetime(startYear,6,21),    datetime(startYear,9,20); % Summer
    datetime(startYear,9,21),    datetime(startYear,12,20) % Autumn
};

%% Loop through each year and season to assign season labels

for thisYear = startYear:endYear
    for iSeason = 1:length(seasons)
        seasonName = seasons{iSeason};

        % Define the start and end dates for the current season and year
        if ~strcmp(seasonName, 'Winter')
            seasonStart = datetime(thisYear, seasonStartEnd{iSeason,1}.Month, seasonStartEnd{iSeason,1}.Day);
        else
            % For winter, the start date is in the previous year
            seasonStart = datetime(thisYear-1, seasonStartEnd{iSeason,1}.Month, seasonStartEnd{iSeason,1}.Day);
        end
        % seasonEnd is adjusted to include the entire last day by adding 
        % one day minus one second to ensure the date range is inclusive of 
        % the end date
        seasonEnd = datetime(thisYear, seasonStartEnd{iSeason,2}.Month, seasonStartEnd{iSeason,2}.Day) + caldays(1) - seconds(1); % Include the last day entirely

        % Assign season labels based on date range
        dataTable.season(dataTable.(dateTimeColName) >= seasonStart &...
            dataTable.(dateTimeColName) <= seasonEnd) = categorical({seasonName});

    end
end

end % addSeason
