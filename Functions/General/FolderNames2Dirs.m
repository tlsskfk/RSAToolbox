function [AllDirNames] = FolderNames2Dirs(DirNames,FileNames)
%Written by David Rothlein
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    FileNames=[];
else
    DirNames=[DirNames,{FileNames}];
end
NumDirNames=size(DirNames,2);
for i = 2:NumDirNames
    TempDirs=[];
    DirNames1=DirNames{1,1};
    if ~iscell(DirNames1)
        DirNames1={DirNames1};
    end
    DirNames2=DirNames{1,i};
    if ~iscell(DirNames2)
        DirNames2={DirNames2};
    end 
    DirNames1=DirNames1(:);
    DirNames2=DirNames2(:);
    for j = 1:length(DirNames1)
        DirName1=DirNames1{j,1};
        for k = 1:length(DirNames2)
            DirName2=DirNames2{k,1};
            if i==NumDirNames && isempty(FileNames)
                TempDirs=[TempDirs;{[DirName1,'/',DirName2,'/']}];
            else
                TempDirs=[TempDirs;{[DirName1,'/',DirName2]}];
            end
        end
    end
    DirNames{1,1}=TempDirs;
end
AllDirNames=DirNames{1,1};        
AllDirNames=strrepCell(AllDirNames,'\','/');
AllDirNames=strrepCell(AllDirNames,'//','/');
end

