function [sceneAverageYearlyMean,sceneAverageMonthlyMean] =...
    calcSceneTimeAveragedData(sceneAverage,time)

%% Calculate monthly and annual average

TT = timetable(time,sceneAverage);
sceneAverageYearlyMean = retime(TT,'yearly','mean');
sceneAverageMonthlyMean = retime(TT,'monthly','mean');

end