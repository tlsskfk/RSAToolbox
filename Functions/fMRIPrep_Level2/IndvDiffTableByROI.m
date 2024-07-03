function [AllROI] = IndvDiffTableByROI(ParcelNames,featureNames,IndvDiffTables,BehTable)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
AllROI = struct;
for j = 1: length(ParcelNames)
ParcelName=ParcelNames{j,1};
ROINames=IndvDiffTables.stat{j,1}.Properties.VariableNames;
ROINames=ROINames(:);
for k = 1:length(ROINames)
ROIName=ROINames{k,1};
dataMat=zeros(size(BehTable,2),length(featureNames));
for i = 1:length(featureNames)
dataMat(:,i)=table2array(IndvDiffTables.stat{j,i}(:,k));
end
dataMat=dataMat';
dataMat=array2table(dataMat,'VariableNames',BehTable.Properties.VariableNames,'RowNames',featureNames);
try
    AllROI.(ParcelName).(ROIName)=dataMat;
catch
    ROIName=strrep(ROIName,'_','');
end
end
end
end