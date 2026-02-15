% data = thisScene;
% figure()
% imagesc(data)
function mask = customAnomalyFilter(data, neighborhoodSize, threshold)
    % Initialize the output mask
    mask = false(size(data));
    
    % Pad the data for boundary cases. A circular padding will avoid
    % repeating anomalously high values at the border
    paddedData = padarray(data, [neighborhoodSize, neighborhoodSize], 'circular');
    
    % Iterate through each pixel
    for i = 1:size(data, 1)
        for j = 1:size(data, 2)
            
            % Extract the local neighborhood
            neighborhood = paddedData(i:i+2*neighborhoodSize, j:j+2*neighborhoodSize);
            
            % Set to 0 values that are NaN (this weill help to remove
            % coloured pixels surrounded by cloud pixels)
            neighborhood(isnan(neighborhood)) = 0;
            
            % Calculate the threshold based on the local neighborhood
            localThreshold = prctile(neighborhood(:), 90); % Adjust percentile as needed
            
            % Check if the central pixel is anomalously high
%             mask(i, j) = data(i, j) > threshold && data(i, j) > localThreshold;
            mask(i, j) = data(i, j) > localThreshold && data(i, j) > threshold;
            
        end
    end
end