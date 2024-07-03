function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = GLM_ComputeTs(fmriprep_table,ExperimentsDir,varargin)
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
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Run AM2 regression for timecourses (default is 1 for yes)
[AM2] = VariableSetter('AM2',1,varargin);
%Run GLM Contrast (default is [] triggering a UI selection)
[RunContrast] = VariableSetter('RunContrast',[],varargin);
%Set batch size for glm. set for optimal tradeoff between speed and memory.
[BatchSize] = VariableSetter('BatchSize',1000,varargin);
% Set Parcellation to run analysis on
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
[LoadActivationPatterns] = VariableSetter('LoadActivationPatterns',[],varargin);
[ActivationPatternName] = VariableSetter('ActivationPatternName',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
[UseAfni] = VariableSetter('UseAfni',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
[BaseTCName] = VariableSetter('SubSamplesName',[],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
Use3DReml = VariableSetter('Use3DReml',0,varargin);
[UseDefaultName] = VariableSetter('UseDefaultName',0,varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_glm,~,glmName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','glm','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_glm.Properties.VariableNames;
        catch
            filePaths_glm=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_glm)
            for j = 1:height(filePaths_glm)
                try
                    load(filePaths_glm.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
                catch
                    continue
                end
                if ~isempty(AnalysisParameters)
                    break
                end
            end
        end
    end
end

if ~isempty(AnalysisParameters)
    ParamNames=fieldnames(AnalysisParameters);
    for i = 1:length(ParamNames)
        if isempty(AnalysisParameters.(ParamNames{i,1}))
            AnalysisParameters.(ParamNames{i,1})='NA_EMPTY';
        end
    end    
    if NoParamEdit==0
        SingleSelect=0; %Allows only a single value to be selected.
        [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect);    
        if ~isempty(EditParams)
            for i = 1:length(EditParams)
                AnalysisParameters.(EditParams{i,1})=[];
            end
        end
    end
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    UseAfni=AnalysisParameters.UseAfni;
    Use3DReml=AnalysisParameters.Use3DReml;
    AfniWorkDir=AnalysisParameters.AfniWorkDir;    
    gmMask=AnalysisParameters.gmMask;
    LoadActivationPatterns=AnalysisParameters.LoadActivationPatterns;
    ConfoundNamesByRun=AnalysisParameters.ConfoundNamesByRun;
    ConfoundNamesBySubj=AnalysisParameters.ConfoundNamesBySubj;
    try
        fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    catch
        fmriprep_table_name=[];
    end
    EventNames=AnalysisParameters.EventNames;
    ActivationPatternName=AnalysisParameters.ActivationPatternName;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    AM2=AnalysisParameters.AM2;
    ConditionNames=AnalysisParameters.ConditionNames;
    ParcelType=AnalysisParameters.ParcelType; 
    ParcelNames=AnalysisParameters.ParcelNames;
    BaseTCName=AnalysisParameters.BaseTCName;    
else
    AnalysisParameters=struct;
    ConditionNames=[];
    ConfoundNamesByRun=[];
    ConfoundNamesBySubj=[];
    EventNames=[];
    AM2Event_Names=[];
    TimeCourseNames=[];
end
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
AnalysisParameters.fmriprep_table_name=fmriprep_table_name;
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants
UseStat='T';
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
if isempty(LoadActivationPatterns)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadActivationPatterns] = uiNameSelect({'Yes','No'},'Load activation patterns?:',SingleSelect);
    if strcmpi(LoadActivationPatterns,'Yes')
        LoadActivationPatterns=1;
    else
        LoadActivationPatterns=0;
    end
end
AnalysisParameters.LoadActivationPatterns=LoadActivationPatterns;
if LoadActivationPatterns==0 
    if isempty(UseAfni)
        SingleSelect=1; %Allows only a single value to be selected.
        [UseAfni] = uiNameSelect({'Yes','No','Hybrid'},'Use AFNI for GLM:',SingleSelect);    
        if strcmpi(UseAfni,'Yes')
            UseAfni=1;
            SingleSelect=1; %Allows only a single value to be selected.
            [Use3DReml] = uiNameSelect({'Yes','No'},'Use AFNI 3dReml:',SingleSelect);
            if strcmpi(Use3DReml,'Yes')
                Use3DReml=1;
            else
                Use3DReml=0;
            end
        elseif strcmpi(UseAfni,'Hybrid')
            UseAfni=2;
        elseif strcmpi(UseAfni,'Hybrid_NoConfounds')
            UseAfni=3;    
        elseif strcmpi(UseAfni,'No')
            UseAfni=0;
        end
    end

    if UseAfni==1 
        if Use3DReml==1
            afniSuffix='_AFNIReml';
        else
            afniSuffix='_AFNI';
        end
        if isempty(AfniWorkDir)
            AfniWorkDir=uigetdir('/','Select AFNI work directory');
            AfniWorkDir=[AfniWorkDir,'/'];
            AfniWorkDir=strrep(AfniWorkDir,'\','/');
        end
    elseif UseAfni==0
        afniSuffix='';
    else
        if isempty(AfniWorkDir)
            AfniWorkDir=uigetdir('/','Select AFNI work directory');
            AfniWorkDir=[AfniWorkDir,'/'];
            AfniWorkDir=strrep(AfniWorkDir,'\','/');
        end    
        afniSuffix='_Hybrid';
    end
end
AnalysisParameters.UseAfni=UseAfni;
AnalysisParameters.Use3DReml=Use3DReml;
AnalysisParameters.AfniWorkDir=AfniWorkDir; 

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    if length(unique(fmriprep_table.session))>1
        useIndicies=find(fmriprep_table.numRuns_bySes)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;        
    else
        useIndicies=find(fmriprep_table.numRuns_bySub)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    end
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='glm';
%Allows you to set name for this particular analysis
if LoadActivationPatterns==0
    %% Compile filepaths for input files for the analysis
    [filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
    filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
    [filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
    filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
    [filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
    filePaths_beh_events=filePaths_beh_events.(beh_eventsName);
end
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
    end
end
if strcmpi(gmMask,'NA_EMPTY')
    gmMask=[];
end
AnalysisParameters.gmMask=gmMask;
if ~isempty(gmMask)
     gmAffix=['_gm',num2str4filename(gmMask,2)];
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end

%% Select variables to include in GLM analysis
% Select confound regressors to include.
if LoadActivationPatterns==0
    ConfoundNames=[];
    if isempty(ConfoundNamesByRun) || (bySS == 1 && isempty(ConfoundNamesBySubj))
        if any(cellfun(@isempty,filePaths_ConfoundTCs)==0)
            for aNum=1:size(filePaths_ConfoundTCs,1)
                if ~isempty(filePaths_ConfoundTCs{aNum,1})            
                    load(filePaths_ConfoundTCs{aNum,1},'ConfoundTCs');            
                    ConfoundNames=unique([ConfoundNames;ConfoundTCs.Properties.VariableNames(:)]);
                end    
            end 
            ConfoundNamesByRun = uiNameSelect(ConfoundNames,'Select within-run confounds to include:');
            if bySS == 1 && UseAfni==0
                ConfoundNamesBySubj = uiNameSelect(ConfoundNames(ismember(ConfoundNames,ConfoundNamesByRun)==0),'Select across-run confounds to include:');
            else
                ConfoundNamesBySubj=[];
            end
        end
    end
    if strcmpi(ConfoundNamesBySubj,'NA_EMPTY')
        ConfoundNamesBySubj=[];
    end
    if strcmpi(ConfoundNamesByRun,'NA_EMPTY')
        ConfoundNamesByRun=[];
    end
end
AnalysisParameters.ConfoundNamesBySubj=ConfoundNamesBySubj;
AnalysisParameters.ConfoundNamesByRun=ConfoundNamesByRun;
%Select Events and Timecourses
EventNamesAll=[];
if LoadActivationPatterns==0
    if any(cellfun(@isempty,filePaths_beh_events)==0)
        if isempty(EventNames) || isempty(TimeCourseNames) || isempty(AM2Event_Names)
            for aNum=1:size(filePaths_beh_events,1)
                if ~isempty(filePaths_beh_events{aNum,1})            
                    load(filePaths_beh_events{aNum,1},'beh_events');            
                    EventNamesAll=unique([EventNamesAll;beh_events.Properties.VariableNames(:)]);
                end    
            end
        end
        if isempty(EventNames)
            EventNames=uiNameSelect(EventNamesAll,'Select events to include');
        elseif strcmpi(EventNames,'NA_EMPTY')
            EventNames=[];
        end
        if isempty(TimeCourseNames)
            TimeCourseNames=uiNameSelect(EventNamesAll,'Select timecourses to include');
        elseif strcmpi(TimeCourseNames,'NA_EMPTY')
            TimeCourseNames=[];
        end        
        if any(AM2==1)  && any(~isempty(TimeCourseNames)) && any(~strcmpi(AM2Event_Names,'NA_EMPTY'))
            for TCNum=1:length(TimeCourseNames)
                [AM2Event_Name] = uiNameSelect([EventNamesAll],['Select AM2 Event for ',TimeCourseNames{TCNum,1},':']);
                AM2Event_Names=[AM2Event_Names;AM2Event_Name];
            end    
        else
            AM2Event_Names=[];
        end
    else
        AM2=0;
    end
end
AnalysisParameters.AM2=AM2;
AnalysisParameters.EventNames=EventNames;
AnalysisParameters.TimeCourseNames=TimeCourseNames;
AnalysisParameters.AM2Event_Names=AM2Event_Names;

%% Set contrast parameters if desired.
if LoadActivationPatterns==0
    if isempty(RunContrast)
        RunContrast=uiEnterName('0',['Run contrast?',newline,'1 = yes, 0 = no']);
        RunContrast=str2num(RunContrast);
    end

    if RunContrast==1
        [cMat,cNames] = uiSelectContrast([EventNames(:);TimeCourseNames(:)]);
    else
        cMat=[];
        cNames=[];
    end

    AnalysisParameters.RunContrast=RunContrast;
    if isempty(ConditionNames)
        ConditionNames = uiNameSelect([EventNames(:);TimeCourseNames(:);cNames(:);unique(AM2Event_Names(:));ConfoundNamesByRun(:);ConfoundNamesBySubj(:)],'Select conditions/contrasts to save:');
    elseif strcmpi(ConditionNames,'NA_EMPTY')
        ConditionNames=[];
    end

    AnalysisParameters.ConditionNames=ConditionNames;
else
    cMat=[];
    cNames=[];
    ConditionNames=[];
end    
if LoadActivationPatterns==1
    [filePaths_ActivationPattern,~,ActivationPatternName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ActivationPatterns','AnalysisName',ActivationPatternName,'TitleTextName','Select activation patterns for RSA:');
    filePaths_ActivationPattern=filePaths_ActivationPattern.WholeBrain;
end
AnalysisParameters.ActivationPatternName=ActivationPatternName;
%Allows you to set name for this particular analysis
if LoadActivationPatterns==1
    UseActivationPatternName=strrep(ActivationPatternName,'ActivationPattern_','GLM_');
    UseActivationPatternName=strrep(UseActivationPatternName,'ActivationPatterns_','GLM_');
    tempName=[UseActivationPatternName,gmAffix];
else
    tempName=strrep(['GLM_',gmAffix,afniSuffix],'__','_');
end
if isempty(AnalysisName)
    tempName=strrep(tempName,'_cl_','_');
    tempName=strrep(tempName,'__','_');
end

if UseDefaultName == 0    
    AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
else
    AnalysisName=tempName;
end
%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%% Parcellation processing
%If necessary, select parcellations to run analyses on. Output of this
%process are the variables:
    %ParcelNames: (N by 1 cell) containing the names of the parcellations 
    %numParcels: Number of parcellations (N) in ParcelNames. 
if ischar(ParcelNames) && ~strcmpi(ParcelNames,'NA_EMPTY')
    ParcelNames={ParcelNames};
end
if isempty(ParcelType)
    ParcelType=uiNameSelect([{'ParcelCoords','GroupParcel','IndvParcel'}],'Select parcellations to run.');
end
AnalysisParameters.ParcelType=ParcelType; 
if strcmpi(ParcelType,'GroupParcel')
    if isempty(ParcelNames)
        ParcelNames=cell(1);
        ParcelDirInfo=dir('Parcellations/');
        count = 1;
        for i = 1:size(ParcelDirInfo,1)
            if ParcelDirInfo(i).isdir==0
                ParcelNames{count,1}=strrep(ParcelDirInfo(i).name,'.mat','');
                count=count+1;
            end
        end
        AllParcelNames=ParcelNames;
        ParcelNames=uiNameSelect([{'WholeBrain'};ParcelNames],'Select parcellations to run.');
    elseif iscell(ParcelNames)
        ParcelNames=ParcelNames(:);
    end
    filePaths_Parcels=[];
elseif strcmpi(ParcelType,'ParcelCoords')
    [ParcelNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','CoordParcels','GetAnalysisNames',1);
    ParcelNames=uiNameSelect([ParcelNames],'Select parcellations to run.');
    AllParcelTable=[];
    for i = 1:length(ParcelNames)
        [tempParcelFileNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','CoordParcels','AnalysisName',ParcelNames{i,1});
        AllParcelTable=[AllParcelTable,tempParcelFileNames];
    end
    ParcelNames=AllParcelTable.Properties.VariableNames(:);
elseif strcmpi(ParcelType,'IndvParcel')
    [AllParcelTable] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','IndvParcels');
    ParcelNames=uiNameSelect([AllParcelTable.Properties.VariableNames],'Select parcellations to run.');
    AllParcelTable=AllParcelTable(:,ParcelNames);
end
   
AnalysisParameters.ParcelNames=ParcelNames;
numParcels=length(ParcelNames);
Parcels=cell(numParcels,1);

%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];
if strcmpi(ParcelType,'GroupParcel')
    if sum(ismember(ParcelNames,'WholeBrain'))==0
        for i = 1:numParcels
            load(['Parcellations/',ParcelNames{i,1}],'UseMask');
            if i ==1            
                analysisMask=single(UseMask>0);
            else
                analysisMask=analysisMask+single(UseMask>0);
            end
            analysisMask=single(analysisMask>0);
        end
    else
       analysisMask=[];
    end
end
BaseAnalysisMask=analysisMask;
if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end

AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.BaseTCName=BaseTCName;
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    analysisMask=BaseAnalysisMask;
    if fmriprep_table.numRuns_bySub(dataInd,1)>0
        subInd=dataInd;
    end
    if fmriprep_table.numRuns_bySes(dataInd,1)>0
        sesInd=dataInd;
    end

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
    descript1='desc-tvals'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(numParcels,1);
    RunParcel=ones(numParcels,1);     
    for parcelNum=1:numParcels
        SavePrefix=[ExperimentsDir,SaveDir,ParcelNames{parcelNum,1},'/'];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
        else  
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
            if exist(SaveNames{parcelNum,1},'file')~=0 && Overwrite==0
                RunParcel(parcelNum,1)=0; 
            end
        end
    end
    if sum(RunParcel(:))==0
        disp(['Skipping-- all files exist: ',SaveName]);
        continue
    end    

    %% Determine number of runs
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    %initialize variable names
    BaseTC=cell(1);
    brain_mask=cell(1);
    RunDur=cell(1);
    Use_ConfoundTCsByRun=cell(1);
    Use_ConfoundTCsBySubj=cell(1);
    Use_Events=cell(1);
    TrialOnsetTimes=cell(1);
    TrialNum=cell(1);
    Use_TimeCourses=cell(1);
    Use_AM2Events=cell(1);
    Use_gmMask=[];    
    TR=fmriprep_table.TR(dataInd,1);
    if LoadActivationPatterns==1
        LoadPath_ActivationPatterns=filePaths_ActivationPattern{dataInd,1};
        if strcmpi(UseStat,'T')
            try
                Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'tVals','brain_mask','DesignMatrix','RegressorNames','AnalysisParameters');
                Use_ActivationPatterns.ActVals=Use_ActivationPatterns.tVals;
                Use_ActivationPatterns.tVals=[];
                brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                continue 
            end
        elseif strcmpi(UseStat,'B')
            try
                Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'bVals','brain_mask','DesignMatrix','RegressorNames','AnalysisParameters');
                Use_ActivationPatterns.ActVals=Use_ActivationPatterns.bVals;
                Use_ActivationPatterns.bVals=[];
                brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                continue 
            end 
        end
        brain_mask=repmat(brain_mask,[numRuns,1]);
    end       
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
        if ~isempty(fmriprep_table.Error{loadInd,1})
            continue
        end    
        %% set load paths and variable names
        %pull load paths
        if LoadActivationPatterns==0 
            LoadPath_Events=filePaths_beh_events{loadInd,1};
            LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
            LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
            %Pull base timecourse data. If it doesn't exist, skip run.
            try
                TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                BaseTC{count,1}=TempLoadData.boldTCs; 
                brain_mask{count,1}=TempLoadData.brain_mask;
                RunDur{count,1}=size(BaseTC{count,1},1);
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
                continue
            end 

            % pull confound TCs
            try
                TempLoadData = load(LoadPath_ConfoundTCs,'ConfoundTCs');
                if ~isempty(ConfoundNamesByRun)
                    ConfoundTCs=TempLoadData.ConfoundTCs(:,TempLoadData.ConfoundTCs.Properties.VariableNames(ismember(TempLoadData.ConfoundTCs.Properties.VariableNames,ConfoundNamesByRun)));
                    Use_ConfoundTCsByRun{count,1}=table2array(ConfoundTCs);
                    Use_ConfoundTCsByRun{count,2}=ConfoundTCs.Properties.VariableNames(:);
                end
                if ~isempty(ConfoundNamesBySubj)
                    ConfoundTCs=TempLoadData.ConfoundTCs(:,TempLoadData.ConfoundTCs.Properties.VariableNames(ismember(TempLoadData.ConfoundTCs.Properties.VariableNames,ConfoundNamesBySubj)));
                    Use_ConfoundTCsBySubj{count,1}=table2array(ConfoundTCs);
                    Use_ConfoundTCsBySubj{count,2}=ConfoundTCs.Properties.VariableNames(:);
                end            
            catch
                disp(['No confound regressors-- input file or variable doesnt exist: ',LoadPath_ConfoundTCs]);
            end   

            %Pull experiment- or  behavior-based events, timecourses and AM2 events
            if ~isempty(LoadPath_Events)
                try
                    TempLoadData = load(LoadPath_Events,'beh_events');
                    if ~isempty(EventNames)
                        Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,EventNames)));
                        Use_Events{count,1}=table2array(Events);
                        Use_Events{count,2}=Events.Properties.VariableNames(:);
                    end 
                    if ~isempty(TimeCourseNames)
                        TimeCourses=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,TimeCourseNames)));
                        Use_TimeCourses{count,1}=table2array(TimeCourses);
                        Use_TimeCourses{count,2}=TimeCourses.Properties.VariableNames(:);
                    end 
                    if ~isempty(AM2Event_Names)
                        AM2Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,AM2Event_Names)));
                        Use_AM2Events{count,1}=table2array(AM2Events);
                        Use_AM2Events{count,2}=AM2Events.Properties.VariableNames(:);
                    end                
                catch
                    disp(['No events-- input file or variable doesnt exist: ',LoadPath_Events]);
                end  
                 try
                    TempLoadData = load(LoadPath_Events,'beh_events');     
                    TrialNum{count,1}=TempLoadData.beh_events.TrialNum; 
                    TrialOnsetTimes{count,1}=TempLoadData.beh_events.TrialOnsetTime; 
                catch
                   disp(['Event load error-- skipping',LoadPath_Events]); 
                   continue
                end           
            end
        end
        if ~isempty(gmMask) && isempty(Use_gmMask)
            LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
            try
                TempLoadData = load(LoadPath_gmMask,'GM_probseg');
                Use_gmMask=single(TempLoadData.GM_probseg > gmMask);
            catch
                disp(['Gray matter mask error-- mask not applied',LoadPath_gmMask]);
                Use_gmMask=(brain_mask{count,1}*0)+1;
            end
        elseif isempty(Use_gmMask)
            Use_gmMask=(brain_mask{count,1}*0)+1;
        end
        
        count=count+1;
    end  
   
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    if isempty(analysisMask)
        analysisMask=Use_gmMask;
    else
        analysisMask=analysisMask.*Use_gmMask;
    end
    %% Run analysis here!!
    try
        if LoadActivationPatterns==0
            [Condition_tVals,brain_mask,RegressorNames,~,DesignMatrix] = GetActivationPatterns(...
                Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
                BaseTC,brain_mask,ConditionNames,...
                'UseAfni',UseAfni,...
                'TrialOnsetTimes',TrialOnsetTimes,...
                'AfniWorkDir',AfniWorkDir,...
                'Use3DReml',Use3DReml,...
                'TR',TR,...
                'TrialNum',TrialNum,...
                'parGLM',1,...
                'ContrastNames',cNames,...
                'Normalize','zscore',...
                'ConfoundTCsBySubj',Use_ConfoundTCsBySubj,...
                'ResampleSizes',RunDur,...
                'AM2',AM2,...
                'AM2_Events',Use_AM2Events,...
                'TimecourseShift',0,...
                'ParcellationMask',analysisMask,...
                'ResampleSizesBoldTC',[],...
                'NormWithinRun',[],...
                'NormAcrossRun','zscore',...
                'Contrasts',cMat,...
                'BatchSize',BatchSize,...
                'ResampleToFit','Y');  
        else
             Condition_tVals=Use_ActivationPatterns.ActVals;
             brain_mask=Use_ActivationPatterns.brain_mask;
             RegressorNames=Use_ActivationPatterns.RegressorNames;
             DesignMatrix=Use_ActivationPatterns.DesignMatrix;
             ConditionNames=Use_ActivationPatterns.AnalysisParameters.ConditionNames;
        end
    catch
        disp('Error')
        continue
    end
    [brainMap] = brainMask2brainMap(Condition_tVals,brain_mask);
    brainSize=size(brainMap);
    for i = 1:numParcels
        
        if strcmpi(ParcelNames{i,1},'WholeBrain')
            brainVals=Condition_tVals;
            save(SaveNames{i,1},'brainVals','brainMap','ConditionNames','brain_mask','RegressorNames','DesignMatrix','AnalysisParameters');
        elseif ~strcmpi(ParcelType,'GroupParcel') 
            tempParcelPaths=table2cell(AllParcelTable(subInd,:));
            load(tempParcelPaths{1,i},'UseMask','UseLabels');
            if unique(UseMask(:)) == 0
                disp(['Skipping ',tempParcelPaths{1,i},' -- ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                continue
            end   
            if ~ismember(brainSize(1,[1:3]),size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize(1,[1:3]),'nearest');
            end    
            UseLabels=UseLabels(:);
            numLabels=length(UseLabels);
            if length(brainSize) ==3
                numCond=1;
            else
                numCond=brainSize(1,4);
            end
            brainVals=nan(numLabels,numCond);
            tempMat=brainMap;
            tempMat=reshape(tempMat,[brainSize(1,1)*brainSize(1,2)*brainSize(1,3),numCond]); 
            UseMask=UseMask(:); 
            for ROI=1:numLabels
                try
                    brainVals(ROI,:)=nanmean(tempMat(UseMask==ROI,:),1);
                catch
                    continue
                end
            end
            try
            brainVals=array2table(brainVals','VariableNames',UseLabels,'RowNames',ConditionNames);
            save(SaveNames{i,1},'brainVals','RegressorNames','DesignMatrix','AnalysisParameters');
            catch
                disp('Error')
                continue
            end            
        else
            load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
            if ~ismember(brainSize(1,[1:3]),size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize(1,[1:3]),'nearest');
            end    
            UseLabels=UseLabels(:);
            numLabels=length(UseLabels);
            if length(brainSize) ==3
                numCond=1;
            else
                numCond=brainSize(1,4);
            end
            brainVals=nan(numLabels,numCond);
            tempMat=brainMap;
            tempMat=reshape(tempMat,[brainSize(1,1)*brainSize(1,2)*brainSize(1,3),numCond]); 
            UseMask=UseMask(:); 
            for ROI=1:numLabels
                try
                    brainVals(ROI,:)=nanmean(tempMat(UseMask==ROI,:),1);
                catch
                    continue
                end
            end
            try
            brainVals=array2table(brainVals','VariableNames',UseLabels,'RowNames',ConditionNames);
            save(SaveNames{i,1},'brainVals','RegressorNames','DesignMatrix','AnalysisParameters');
            catch
                disp('Error')
                continue
            end
        end
    end        
    toc
end
ConditionNames=ConditionNames(:);
for parcelNum=1:numParcels  
    ParcelName=ParcelNames{parcelNum,1};
    parcelName=ParcelName;
    if strcmpi(ParcelName,'WholeBrain')
        continue
    end        
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
    end
    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);     
    end
    if ~exist(GroupFiguresDir,'file')
        mkdir(GroupFiguresDir);     
    end
    if ~exist(GroupBrainMapsDir,'file')
        mkdir(GroupBrainMapsDir);     
    end     
    try
        load(['Parcellations/',parcelName],'UseLabels','UseMask');
        UseLabels=UseLabels(:);
    catch
        try
            [compiledData] = CompileND(ExperimentsDir,fmriprep_table,'AnalysisType','CoordParcels','LoadVarName','ROInames','LoadVarFormat','Cell','DataType','Other','AnalysisName',parcelName,'SubjectOrRun',SubjectOrRun,'TableOrArray','Table','VertOrMat','Matrix','ColumnCompile','All');
            UseLabels=compiledData(:,1,1);
            UseMask=[];
        catch
            load(['Parcellations/IndvParcels/',fmriprep_table_name,'/',parcelName],'UseLabels','UseMask')
        end
    end
    save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
    if isempty(ConditionNames)
        [filePaths_ConditionNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType',AnalysisType,'AnalysisName',AnalysisName);
        filePaths_ConditionNames=filePaths_ConditionNames.(ParcelName);  
        for j = 1:length(filePaths_ConditionNames)
            try
                load(filePaths_ConditionNames{j,1},'brainVals');
                ConditionNames=brainVals.Properties.RowNames;
            catch
                continue
            end
            if ~isempty(ConditionNames)
                break
            end
        end
    end
    for condNum= 1:length(ConditionNames)
        condName=ConditionNames{condNum,1};
        [vecVals] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','brainVals',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Array',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName,...
            'RowCompile',condName);
        [GroupActivationResults,ReliabilityResults] = GroupVecParcellationSummaryFigures(vecVals,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,'SavePrefix',condName); 
        Activation_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
        Activation_losoReliability=ReliabilityResults.Vec_losoReliability;        
        save([GroupAnalysisDir,'GroupActivationResults_',condName],'GroupActivationResults');
        save([GroupAnalysisDir,'Group_ParcelReliability_',condName],'Activation_SplitHalfReliability','Activation_losoReliability');   
    end
end  
    
end 