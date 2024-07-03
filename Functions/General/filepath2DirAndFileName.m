function [outDirs,fileNames] = filepath2DirAndFileName(inCell)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if istable(inCell)
    outDirs=inCell;
    fileNames=inCell;    
    inCell=table2cell(inCell);
    [outDirs(:,:),fileNames(:,:)]=cellfun(@filepathSplitter,inCell,'UniformOutput',0);
else
    [outDirs,fileNames]=cellfun(@filepathSplitter,inCell,'UniformOutput',0);
end

end

function [outDir,fileName] = filepathSplitter(inName)
if isempty(inName)
    outDir=[];
    fileName=[];
else
    if iscell(inName)
        inName=inName{1,1};
    end
    splitPt=strfind(inName,'/');
    splitPt=splitPt(end);
    outDir=inName(1,1:splitPt);
    fileName=inName(1,splitPt+1:end);
end
end