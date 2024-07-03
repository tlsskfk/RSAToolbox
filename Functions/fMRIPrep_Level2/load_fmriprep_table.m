function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table(varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ExperimentsDir] = VariableSetter('ExperimentsDir',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);

if isempty(ExperimentsDir)
    ExperimentsDir=uigetdir('/','Select Experiments Directory');
    ExperimentsDir=[ExperimentsDir,'/'];
end
ExperimentsDir=strrep(ExperimentsDir,'\','/');
tempDir=strrep(ExperimentsDir,'/Experiments/','/');
slashInd=strfind(tempDir,'/');
ExpName=tempDir(1,slashInd(1,end-1):slashInd(1,end));
if isempty(fmriprep_table_name)
    [~,fmriprep_table_names]=getFolderAndFileNames(['fmriprep_table',ExpName]);
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect([{'none- make table'};fmriprep_table_names],'Select fmriprep_table used:',SingleSelect);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
end

if strcmpi(fmriprep_table_name,'none- make table')
    fmriprep_table=[];
else
    try
        load(['fmriprep_table',ExpName,fmriprep_table_name]);
    catch
        fmriprep_table=[];
    end
end
end
