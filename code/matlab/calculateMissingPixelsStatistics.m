function [fracMissingPixels,isValidImage,fracMissingPixelsByMonth,...
    imageDistribByMonthAndValidity,M] = calculateMissingPixelsStatistics(...
    chlaSatellite,timeVector,nImages,nPixelsInImage,fracMissingPixelToConsiderInvalidImage)

% Calculate the percentage of images where > 98% pixels are missing
% (invalid images)
fracMissingPixels = NaN(1,nImages);
isValidImage = NaN(1,nImages);
for iImage = 1:nImages
    chlaArray = squeeze(chlaSatellite(:,:,iImage)); 
    chlaColumn = reshape(chlaArray, [], 1);
    nNaNpixels = sum(isnan(chlaColumn));
    fracMissingPixels(iImage) = nNaNpixels/nPixelsInImage;
    if (fracMissingPixels(iImage) > fracMissingPixelToConsiderInvalidImage/100)
        isValidImage(iImage) = 0;
    else
        isValidImage(iImage) = 1;
    end
end

% Grouping by month and calculating frequencies of valid vs invalid
M = table(timeVector,isValidImage',fracMissingPixels',...
    'VariableNames',{'date','valid','fracMissingPixels'});
M.month = month(M.date);
M.year = year(M.date);
M.season = strings(size(M,1),1); 

imageDistribByMonthAndValidity = zeros(12,2);
fracMissingPixelsByMonth = zeros(12,1);
for iMonth = 1:12
    thisMonthData = M.valid(M.month(:) == iMonth);
    imageDistribByMonthAndValidity(iMonth,1) = sum(thisMonthData == 1); % valid images
    imageDistribByMonthAndValidity(iMonth,2) = sum(thisMonthData == 0); % invalid imags
    thisMonthValidData = M.fracMissingPixels(M.month(:) == iMonth & M.valid == 1);
    fracMissingPixelsByMonth(iMonth) = mean(thisMonthValidData,'omitnan');
end

% Assign season to data
seasons = {'Winter','Spring','Summer','Autumn'};
seasonStartEnd = {
    datetime('1997-12-21'), datetime('1998-03-20');
    datetime('1998-03-21'), datetime('1998-06-20');
    datetime('1998-06-21'), datetime('1998-09-20');
    datetime('1998-09-21'), datetime('1998-12-20');
};

for thisYear = min(M.year):max(M.year)
    for iSeason = 1:length(seasons)

        if ~strcmp(seasons{iSeason},'Winter') 
            seasonStart = datetime(thisYear, seasonStartEnd{iSeason,1}.Month, seasonStartEnd{iSeason,1}.Day);
        elseif strcmp(seasons{iSeason},'Winter') % in the winter, start the season in the year before
            seasonStart = datetime(thisYear-1, seasonStartEnd{iSeason,1}.Month, seasonStartEnd{iSeason,1}.Day);
        end
        seasonEnd = datetime(thisYear, seasonStartEnd{iSeason,2}.Month, seasonStartEnd{iSeason,2}.Day);

        % Extract data for the current season and year
        iCurrYearAndSeason = find(M.date >= seasonStart & M.date <= seasonEnd);
        if (~isempty(iCurrYearAndSeason))
            switch seasons{iSeason}
                case 'Winter'
                    M.season(iCurrYearAndSeason) = 'Winter';
                case 'Spring'
                    M.season(iCurrYearAndSeason) = 'Spring';
                case 'Summer'
                    M.season(iCurrYearAndSeason) = 'Summer';
                case 'Autumn'
                    M.season(iCurrYearAndSeason) = 'Autumn';
            end
        end
    end
end
    
end