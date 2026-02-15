function [regriddedData, interpolatedData] = regridAndInterpolateData(dataset,...
    dataLatVec,dataLonVec,dataTimeVec,dataDepthVec,areaStudyLatVec,areaStudyLonVec)

% Determine if the dataset has a depth dimension
hasDepthDimension = ~isempty(dataDepthVec);
    
% Set up query points for regridding based on whether there's a depth dimension
if hasDepthDimension
    % Set up query points for regridding
    [qX, qY, qT, qD] = ndgrid(areaStudyLatVec, areaStudyLonVec, dataTimeVec, dataDepthVec);
    % Set up the original grid points
    [X, Y, T, D] = ndgrid(dataLatVec, dataLonVec, dataTimeVec, dataDepthVec);
    % Create the gridded interpolant using the input dataset
    F = griddedInterpolant(X, Y, T, D, dataset);
    % Evaluate the interpolant on the new grid to obtain regridded data
    regriddedData = F(qX, qY, qT, qD);
else
    % Set up query points for regridding
    [qX, qY, qT] = ndgrid(areaStudyLatVec, areaStudyLonVec, dataTimeVec);
    % Set up the original grid points
    [X, Y, T] = ndgrid(dataLatVec, dataLonVec, dataTimeVec);
    % Create the gridded interpolant using the input dataset
    F = griddedInterpolant(X, Y, T, dataset);
    % Evaluate the interpolant on the new grid to obtain regridded data
    regriddedData = F(qX, qY, qT);
end
    
% Interpolate data to fill in any gaps
interpolatedData = Racault2014interpolation(regriddedData);

end