function [AnalysisParameters] = GLM_ComputeResids(fmriprep_table,ExperimentsDir,varargin)
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
[SubjectOrRun] = VariableSetter('SubjectOrRun',['Run'],varargin);
%Run AM2 regression for timecourses (default is 1 for yes)
[AM2] = VariableSetter('AM2',1,varargin);
%Set batch size for glm. set for optimal tradeoff between speed and memory.
[BatchSize] = VariableSetter('BatchSize',1000,varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[BaseTCName] = VariableSetter('SubSamplesName',[],varargin);
[UseAfni] = VariableSetter('UseAfni',[],varargin);
[Use3DReml] = VariableSetter('Use3DReml',0,varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_resid,~,residName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','func','SubjectOrRun','Run','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_resid.Properties.VariableNames;
        catch
            filePaths_resid=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_resid)
            for j = 1:height(filePaths_resid)
                try
                    load(filePaths_resid.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
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
    SingleSelect=0; %Allows only a single value to be selected.
    [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect);    
    if ~isempty(EditParams)
        for i = 1:length(EditParams)
            AnalysisParameters.(EditParams{i,1})=[];
        end
    end
    SubjectOrRun='Run';
    try
        UseAfni=AnalysisParameters.UseAfni;
        AfniWorkDir=AnalysisParameters.AfniWorkDir;
        Use3DReml=AnalysisParameters.Use3DReml;
    catch
        UseAfni=0;
        AfniWorkDir=[];
        Use3DReml=0;
    end   
    ConfoundNames=AnalysisParameters.ConfoundNames;
    EventNames=AnalysisParameters.EventNames;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    AM2=AnalysisParameters.AM2;
    BaseTCName=AnalysisParameters.BaseTCName;    
else
    AnalysisParameters=struct;
    ConfoundNames=[];
    EventNames=[];
    AM2Event_Names=[];
    TimeCourseNames=[];
end

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 

if isempty(UseAfni)
    SingleSelect=1; %Allows only a single value to be selected.
    [UseAfni] = uiNameSelect({'Yes','No','Hybrid'},'Use AFNI for GLM:',SingleSelect);    
    if strcmpi(UseAfni,'Yes')
        UseAfni=1;
        Use3DReml = 'No';
        if strcmpi(Use3DReml,'Yes')
            Use3DReml=1;
        else
            Use3DReml=0;
        end
    elseif strcmpi(UseAfni,'Hybrid')
        UseAfni=2;
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
AnalysisParameters.UseAfni=UseAfni;
AnalysisParameters.Use3DReml=Use3DReml;
AnalysisParameters.AfniWorkDir=AfniWorkDir; 

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='func';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['ResidTC',afniSuffix,'_',genDateString],['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%% Compile filepaths for input files for the analysis
restingData=0;
try
    [filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
    filePaths_beh_events=filePaths_beh_events.(beh_eventsName);
catch
    disp('No behavior found! Okay if resting data.')
    filePaths_beh_events=[];
    restingData=1;
end
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);  
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
AnalysisParameters.BaseTCName=BaseTCName;

%% Select variables to include in GLM analysis
% Select confound regressors to include.
%% Select variables to include in GLM analysis
% Select confound regressors to include.
if isempty(ConfoundNames)
    if any(cellfun(@isempty,filePaths_ConfoundTCs)==0)
        for aNum=1:size(filePaths_ConfoundTCs,1)
            if ~isempty(filePaths_ConfoundTCs{aNum,1})            
                load(filePaths_ConfoundTCs{aNum,1},'ConfoundTCs');            
                ConfoundNames=unique([ConfoundNames;ConfoundTCs.Properties.VariableNames(:)]);
            end    
        end 
        ConfoundNames = uiNameSelect(ConfoundNames,'Select within-run confounds to include:');
    end
end
if strcmpi(ConfoundNames,'NA_EMPTY')
    ConfoundNames=[];
end
AnalysisParameters.ConfoundNames=ConfoundNames;
%Select Events and Timecourses
EventNamesAll=[];
if restingData==0
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
        if AM2==1  && ~isempty(TimeCourseNames) && ~strcmpi(AM2Event_Names,'NA_EMPTY')
            for TCNum=1:length(TimeCourseNames)
                [AM2Event_Name] = uiNameSelect([EventNamesAll],['Select AM2 Event for ',TimeCourseNames{TCNum,1},':']);
                AM2Event_Names=[AM2Event_Names;AM2Event_Name];
            end    
        else
            AM2Event_Names=[];
        end
    end
else
    AM2=0;
end
AnalysisParameters.AM2=AM2;
AnalysisParameters.EventNames=EventNames;
AnalysisParameters.TimeCourseNames=TimeCourseNames;
AnalysisParameters.AM2Event_Names=AM2Event_Names;

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

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    %% Set save directory and save name
    SaveDir=strrep(fmriprep_table.funcDir{dataInd,1},'/func/',['/',AnalysisType,'/',AnalysisName,'/']);
    SaveDir=strrep(SaveDir,'/fmriprep/','/matlab/');
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-boldTCs'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis does not involves parcellations:
    SavePrefix=[ExperimentsDir,SaveDir,'/'];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
    end    
    SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
        disp(['Skipping-- file exists: ',SaveName]);
        continue
    end    
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
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
        
        if ~isempty(filePaths_beh_events)
            LoadPath_Events=filePaths_beh_events{loadInd,1};
        else
            LoadPath_Events=[];
        end
        LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
        LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
        
        %initialize variable names
        BaseTC=cell(1);
        brain_mask=cell(1);
        RunDur=cell(1);
        Use_ConfoundTCs=cell(1);
        Use_Events=cell(1);
        Use_TimeCourses=cell(1);
        Use_AM2Events=cell(1);
        TrialOnsetTimes=cell(1);
        TrialNum=cell(1);        
        TR=fmriprep_table.TR(dataInd,1);
        %Pull base timecourse data. If it doesn't exist, skip run.
        try
            variableInfo = who('-file', LoadPath_BaseTC);
            if any(contains(variableInfo,'boldTCs'))
                TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                BaseTC{count,1}=TempLoadData.boldTCs;
            elseif any(contains(variableInfo,'Resids'))   
                TempLoadData = load(LoadPath_BaseTC,'Resids','brain_mask');
                BaseTC{count,1}=TempLoadData.Resids;
            end    
            brain_mask{count,1}=TempLoadData.brain_mask;
            RunDur{count,1}=size(BaseTC{count,1},1);
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
            continue
        end         
        
        % pull confound TCs
        try
            TempLoadData = load(LoadPath_ConfoundTCs,'ConfoundTCs');
            if ~isempty(ConfoundNames)
                ConfoundTCs=TempLoadData.ConfoundTCs(:,TempLoadData.ConfoundTCs.Properties.VariableNames(ismember(TempLoadData.ConfoundTCs.Properties.VariableNames,ConfoundNames)));
                Use_ConfoundTCs{count,1}=table2array(ConfoundTCs);
                Use_ConfoundTCs{count,2}=ConfoundTCs.Properties.VariableNames(:);
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
                TrialNum{count,1}=TempLoadData.beh_events.TrialNum; 
                TrialOnsetTimes{count,1}=TempLoadData.beh_events.TrialOnsetTime; 
            catch
                disp(['No events-- input file or variable doesnt exist: ',LoadPath_Events]);
            end          
        end
        count=count+1;
    end  
   
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    try
    %% Run analysis here!!
    [Resids,brain_mask,~,DesignMatrix,AfniInfo] = GetResids(...
        Use_Events,Use_TimeCourses,Use_ConfoundTCs,...
        BaseTC,brain_mask,...
        'UseAfni',UseAfni,...
        'TrialOnsetTimes',TrialOnsetTimes,...
        'AfniWorkDir',AfniWorkDir,...
        'TR',TR,...
        'TrialNum',TrialNum,...
        'parGLM',1,...
        'Normalize','zscore',...
        'ResampleSizes',RunDur,...
        'AM2',AM2,...
        'AM2_Events',Use_AM2Events,...
        'TimecourseShift',0,...
        'ResampleSizesBoldTC',[],...
        'NormWithinRun',[],...
        'NormAcrossRun','zscore',...
        'BatchSize',BatchSize,...
        'ResampleToFit','Y');      
 
    save(SaveNames{1,1},'Resids','brain_mask','DesignMatrix','AfniInfo','AnalysisParameters'); 
    catch
        delete([AfniWorkDir,'*']);        
        disp('Error')
        continue        
    end
    toc
end
end 