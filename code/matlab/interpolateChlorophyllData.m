function interpolatedData = interpolateChlorophyllData(data)

% Linear interpolation of data as in Racault et al. (2014). We have a 3D 
% data array of chlorophyll a, the dimensions are latitude x longitude x time. 
% We want to construct a function that reduces the number of missing data 
% by doing a spatial and temporal linear interpolation. The interpolation 
% scheme has to be applied sequentially in the order: longitude, latitude, 
% and time using a three-point window. If one of the points bordering the 
% gap along the indicated axis is invalid, it will be omitted from the 
% calculation, whilst if two surrounding points are invalid, then the 
% gap will not be filled. 

    [latSize, lonSize, timeSize] = size(data);
    interpolatedData = interpolateLongitude(data,latSize,timeSize);
    interpolatedData = interpolateLatitude(interpolatedData,lonSize,timeSize);
    interpolatedData = interpolateTime(interpolatedData,latSize,lonSize);
end

% Local functions to this script

function interpolatedData = interpolateLongitude(inputData,latSize,timeSize)
    interpolatedData = inputData;
    for iTime = 1:timeSize
        for iLat = 1:latSize
            interpolatedData(iLat, :, iTime) = interpolateDimension(inputData(iLat, :, iTime));
        end
    end
end

function interpolatedData = interpolateLatitude(inputData,lonSize,timeSize)
    interpolatedData = inputData;
    for iTime = 1:timeSize
        for iLon = 1:lonSize
            interpolatedData(:, iLon, iTime) = interpolateDimension(inputData(:, iLon, iTime));
        end
    end
end

function interpolatedData = interpolateTime(inputData,latSize,lonSize)
    interpolatedData = inputData;
    for iLat = 1:latSize
        for iLon = 1:lonSize
            interpolatedData(iLat, iLon, :) = interpolateDimension(inputData(iLat, iLon, :));
        end
    end
end

function interpolatedValues = interpolateDimension(values)
    interpolatedValues = values;
    for i = 2:(length(values) - 1)
        window = values(i - 1:i + 1);
        if isnan(interpolatedValues(i))
            if any(~isnan(window))
                interpolatedValues(i) = mean(window, 'omitnan');
            end
        end
    end
end


% % Dimensions of the input data
% [latSize, lonSize, timeSize] = size(data);
% 
% % Initialize the interpolated data
% interpolatedDataFirstRound = data;
%         
% for iTime = 1:timeSize
%     
%     for iLat = 1:latSize
%         
%         % Iterate over the longitude axis
%         for iLon = 2:(lonSize - 1)
%             % Create a window of three consecutive latitude points
%             window = data(iLat, iLon - 1:iLon + 1, iTime);
% 
%             % Check if the center point is NaN
%             if isnan(interpolatedDataFirstRound(iLat, iLon, iTime))
%                 % Check if at least one of the bordering points is valid
%                 if any(~isnan(window(:)))
%                     % Perform interpolation by taking the mean
%                     interpolatedDataFirstRound(iLat, iLon, iTime) = mean(window, 'omitnan');
%                 end
%             end
%         end
%         
%     end
% end
% 
% interpolatedDataSecondRound = interpolatedDataFirstRound;
% 
% for iLon = 1:lonSize
% 
%     for iTime = 1:timeSize
%     
%         for iLat = 2:(latSize - 1)
%             % Create a window of three consecutive latitude points
%             window = interpolatedDataFirstRound(iLat - 1:iLat + 1, iLon, iTime);
% 
%             % Check if the center point is NaN
%             if isnan(interpolatedDataSecondRound(iLat, iLon, iTime))
%                 % Check if at least one of the bordering points is valid
%                 if any(~isnan(window(:)))
%                     % Perform interpolation by taking the mean
%                     interpolatedDataSecondRound(iLat, iLon, iTime) = mean(window, 'omitnan');
%                 end
%             end
%         end
%         
%     end
% end
%             
% interpolatedDataThirdRound = interpolatedDataSecondRound;            
% 
% for iLat = 1:latSize
%     
%     for iLon = 1:lonSize
% 
%         for iTime = 2:(timeSize - 1)
%     
%             % Create a window of three consecutive latitude points
%             window = interpolatedDataSecondRound(iLat, iLon, iTime - 1:iTime + 1);
% 
%             % Check if the center point is NaN
%             if isnan(interpolatedDataThirdRound(iLat, iLon, iTime))
%                 % Check if at least one of the bordering points is valid
%                 if any(~isnan(window(:)))
%                     % Perform interpolation by taking the mean
%                     interpolatedDataThirdRound(iLat, iLon, iTime) = mean(window, 'omitnan');
%                 end
%             end
%         end
%         
%     end
% end
% 
