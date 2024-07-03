function [OutTable,OutNumTable,OutOthTable,OutCellTable] = AverageTableBySS(DataTable,GroupingVar)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
OutNumTable=[];
OutOthTable=[];
OutCellTable=[];
%Determine data types
tempT=table2cell(DataTable(1,:));
DataTypes=cellfun(@class,tempT,'UniformOutput',false);
NumInd=ismember(DataTypes,{'single','double','int8','int16','int32','int64','uint8','uint16','uint32','uint64'});
numTable=DataTable(:,NumInd);
cellInd=table2array(varfun(@iscell,numTable));
cellTable=numTable(:,cellInd);
numTable(:,cellInd)=[];
numVarsNames=numTable.Properties.VariableNames;
othTable=DataTable(:,NumInd==0);

%Average over numTable and select over DataTable
if ischar(GroupingVar)
    GroupingVar=DataTable.(GroupingVar);
end

GroupNames=unique(GroupingVar);

for i = 1:length(GroupNames)
    selectInd=ismember(GroupingVar,GroupNames{i,1});
    tempNumTable=array2table(nanmean(table2array(numTable(selectInd,:)),1),'VariableNames',numVarsNames);
    OutNumTable=[OutNumTable;tempNumTable];
    tempOthTable=othTable(selectInd,:);
    tempOthTable=tempOthTable(1,:);
    OutOthTable=[OutOthTable;tempOthTable];
    tempCellTable=cellTable(selectInd,:);
    tempCellTable=tempCellTable(1,:);
    OutCellTable=[OutCellTable;tempCellTable];
    
end    
OutTable=[OutOthTable,OutCellTable,OutNumTable];
end

