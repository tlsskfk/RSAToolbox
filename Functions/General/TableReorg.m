function [OutTables,UseParameters] = TableReorg(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[UseDir] = VariableSetter('UseDir',[],varargin);
[TableFileNames] = VariableSetter('TableFileNames',[],varargin);
[FileVarNames] = VariableSetter('FileVarNames',[],varargin);
[TableVarNames] = VariableSetter('TableVarNames',[],varargin);
[OutTableVarNames] = VariableSetter('OutTableVarNames',[],varargin);
[UseParameters] = VariableSetter('UseParameters',[],varargin);

if ~isempty(UseParameters)
    UseDir=UseParameters.UseDir;
    TableFileNames=UseParameters.TableFileNames;
    FileVarNames=UseParameters.FileVarNames;
    TableVarNames=UseParameters.TableVarNames;
    OutTableVarNames=UseParameters.OutTableVarNames;
end

OutTables=struct;

if isempty(UseDir)
    UseDir=uigetdir();
    UseDir=[UseDir,'/'];
end
if isempty(TableFileNames)
    [~,TableFileNames]=getFolderAndFileNames(UseDir);
    singleselect=0;
    [TableFileNames] = uiNameSelect(TableFileNames,'Select Files to Include:',singleselect);
end
if isempty(FileVarNames)||isempty(TableVarNames)

tempFile=load([UseDir,TableFileNames{1,1}])

end