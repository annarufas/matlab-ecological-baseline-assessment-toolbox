function [sceneAverage] = calcScenePonderedDepthAverage(dataset,...
    time,nDepthLevels,nTimeSteps)

sceneAverageByDepth = NaN(nTimeSteps,nDepthLevels);
sceneNumPixelsByDepth = NaN(nTimeSteps,nDepthLevels);
for iDepthLevel = 1:nDepthLevels
    [sceneAverageByDepth(:,iDepthLevel),sceneNumPixelsByDepth(:,iDepthLevel)] =...
        calcSceneAverage(dataset,time);
end

% Weighted depth mean
sceneAverage = NaN(nTimeSteps,1);
for iTimeStep = 1:nTimeSteps
    totNumPixels = sum(sceneNumPixelsByDepth(iTimeStep,:),'omitnan');
    accum = 0;
    for iDepthLevel = 1:nDepthLevels
        if (~isnan(sceneAverageByDepth(iTimeStep,iDepthLevel)))
            accum = accum + sceneAverageByDepth(iTimeStep,iDepthLevel)...
                .*sceneNumPixelsByDepth(iTimeStep,iDepthLevel);
        end
    end
    sceneAverage(iTimeStep) = accum/totNumPixels;
end
   
end % calcScenePonderedDepthAverage