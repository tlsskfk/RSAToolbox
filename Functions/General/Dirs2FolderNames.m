function [FolderNames] = Dirs2FolderNames(PathList)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numPaths=size(PathList,1);

for i = 1:numPaths
    dInd=strfind(PathList{i,1},'/');
    if dInd(end) == length(PathList{i,1})
        dInd(:,end)=[];
        PathList{i,1}(:,end)=[];
    end
    if dInd(1,1) == 1
        dInd(:,1)=[];
        PathList{i,1}(:,1)=[];
    end
    for j= 0:length(dInd)
        if j==0
            snip=PathList{i,1}(1,1:dInd(1,1)-1);
        elseif j==length(dInd)
            snip=PathList{i,1}(1,dInd(1,j)+1:end);
        else
            snip=PathList{i,1}(1,dInd(1,j)+1:dInd(1,j+1)-1);
        end
        FolderNames{i,j+1}=snip;
    end   
end


end

