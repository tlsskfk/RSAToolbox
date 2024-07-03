% Get the structure array containing the file information
fileInfo = dir('Raw_RSM/');
% Extract the names of the files into a cell array
fileNames = {fileInfo.name};
% Convert the cell array to a character array
fileNamesCharArray = char(fileNames);
% Display the character array
disp(fileNamesCharArray);
% Filter out directories
fileNames = {fileInfo(~[fileInfo.isdir]).name};
% Convert to character array
fileNames = string(fileNames);
structNames=erase(fileNames, '.mat');
for i=1:length(fileNames)
       input=load(fullfile('Raw_RSM', fileNames(i)));
       name=string(fieldnames(input));
       input=input.(name);
       [pRSMs] = WordStudy_pRMSConvert(input, 1);
       save(fullfile('pRSMs/', strcat(structNames(i), '_pRSM', '.mat')), "pRSMs")
end

