function [folderNames,fileNames]=getFolderAndFileNames(dirName,IncludeStr)
%Written by David Rothlein
% Complies vertical cell arrays of file and folder names contained 
% within a directory.

% Input: 
%   dirName - a string specifying the name of the directory to search
%   IncludeStr -  a string that will be used to filter output names to files or
%                 folder names that include text in IncludeStr. Default is
%                 [] which returns all file and folder names
% Other details:
% Automatically removes '.' and '..' folder names.
% If no files or folders are found, an empty array is returned.

if nargin == 1
    IncludeStr=[];
end

dirInfo=dir(dirName);
folderNames=cell(1);
fileNames=cell(1);
folderCount=1;
fileCount=1;
if ~isempty(IncludeStr)
    for j = 1:size(dirInfo,1)
        if dirInfo(j).isdir==1 && contains(dirInfo(j).name,IncludeStr)
            folderNames{folderCount,1}=dirInfo(j).name;
            folderCount=folderCount+1;
        elseif dirInfo(j).isdir==0 && contains(dirInfo(j).name,IncludeStr)
            fileNames{fileCount,1}=dirInfo(j).name;
            fileCount=fileCount+1;
        end
    end
else    
    for j = 1:size(dirInfo,1)
        if strcmp(dirInfo(j).name,'.') || strcmp(dirInfo(j).name,'..')
            continue
        end
        if dirInfo(j).isdir==1
            folderNames{folderCount,1}=dirInfo(j).name;
            folderCount=folderCount+1;
        else
            fileNames{fileCount,1}=dirInfo(j).name;
            fileCount=fileCount+1;
        end
    end    
end

if iscell(folderNames)
    if isempty(folderNames{1})
        folderNames=[];
    end
end
if iscell(fileNames)
    if isempty(fileNames{1})
        fileNames=[];
    end
end
% if folderCount==1
%     folderNames=[];
% end
% if fileCount==1
%     fileNames=[];
% end

end
