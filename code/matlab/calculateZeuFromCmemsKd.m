function [cmemsZeu] = calculateZeuFromCmemsKd(cmemsData)

% Notice that the light extinction coefficient is provided for different
% depth levels. Kd is the slope of a curve so it can only be one value.
% Therefore, we will take the first value provided.

% Find idx of kd product in the CMEMS array
idxKdCmems = find(strcmp({cmemsData.ID}, 'mod_bgc_reg_kd'));
cmemsKd = cmemsData(idxKdCmems).dataset;

% Initialise zeu array
% The size of cmemsZeu matches the first three dimensions (lat, lon, time) of cmemsKd
cmemsZeu = NaN(size(cmemsKd,1:3));
for iTimeStep = 1:numel(cmemsData(idxKdCmems).time)
    for iLon = 1:numel(cmemsData(idxKdCmems).lon)
        for iLat = 1:numel(cmemsData(idxKdCmems).lat)
            % Extract the kd value (first value in depth dimension
            kdValue = cmemsKd(iLat,iLon,iTimeStep,1);
            if (kdValue > 0 && ~isnan(kdValue))
                % Use a standard formula to calculate zeu, zeu = 4.6 / kd
                cmemsZeu(iLat,iLon,iTimeStep) = 4.6./kdValue; 
            end
        end
    end
end

end