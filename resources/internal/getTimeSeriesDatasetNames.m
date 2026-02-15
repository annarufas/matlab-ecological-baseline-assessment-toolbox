function [datasetNames] = getTimeSeriesDatasetNames(datasetTable)
    populatedIDs = ~cellfun(@isempty, {datasetTable.ID});
    datasetNames = {datasetTable(populatedIDs).ID};
end