function [AnalysisParameters] = SaveIndv_Coords(fmriprep_table,ExperimentsDir,CoordLabels,Coords,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end
%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[''],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',['Subject'],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

ExperimentsDir=strrep(ExperimentsDir,'\','/');
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
CoordLabels=strrepCell(CoordLabels,'-','_');
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ROICoords';
%Allows you to set name for this particular analysis
ROInames=[];
if isnumeric(Coords)
    Coords=num2cell(Coords,2);
end
for i = 1:height(fmriprep_table)
    subStr=['_',fmriprep_table.sub{i,1},'_'];
    tempLabels=CoordLabels(contains(CoordLabels,subStr),:);
    tempLabels=strrepCell(tempLabels,subStr,'_');
    ROInames=unique([ROInames;tempLabels]);
end
%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=dataInd;
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['ROICoords'],['Enter name for ',AnalysisType,newline,'analysis below:']);
end
AnalysisParameters.ROInames=ROInames;
AnalysisParameters.SubjectOrRun=SubjectOrRun;

iniPercentComplete=0; %Used to display progress
ROIInfo=[];
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    %% Set save directory and save name
    SaveDir=[fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
        SaveName=strrep(SaveName,'_run-1','');
    end
    descript1='desc-ROICoords'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(1);
    SavePrefix=[ExperimentsDir,SaveDir,'/']; 
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    else  
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
            disp(['Skipping-- all files exist: ',SaveName]);
        end
    end
    %% Determine number of runs
    subStr=['_',fmriprep_table.sub{dataInd,1},'_'];
    ssSelectInd = contains(CoordLabels,subStr);
    tempLabels=CoordLabels(ssSelectInd,:);
    tempLabels=strrepCell(tempLabels,subStr,'_');
    tempCoords=Coords(ssSelectInd,:);
    ROICoords=cell2table(tempCoords','VariableNames',tempLabels);
    save(SaveNames{1,1},'ROICoords');      
end
