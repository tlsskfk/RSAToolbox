function [outPaths] = sub2runFilePath(numRuns_bySub,filePaths)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here;
outPaths=cell(size(numRuns_bySub,1),1);
tablePaths=[];
if istable(filePaths)
    tablePaths=filePaths;
    filePaths=table2cell(filePaths);   
end
inds=find(numRuns_bySub)';
for i = inds
    usePath=filePaths{i,1};
    for j = 1:numRuns_bySub(i,1)
        outPaths{i+j-1,1}=usePath;
    end
end
if ~isempty(tablePaths)
    outPaths=cell2table(outPaths,'VariableNames',tablePaths.Properties.VariableNames);
end
end