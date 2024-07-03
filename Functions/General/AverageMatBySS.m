function [OutData] = AverageMatBySS(Data,GroupingVar,Dim)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    Dim = ndims(Data);
end
GroupNames=unique(GroupingVar);
OutData=Data*0;
if Dim == 1

for i = 1:length(GroupNames)
    selectInd=ismember(GroupingVar,GroupNames{i,1});
    tempData=Data
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

