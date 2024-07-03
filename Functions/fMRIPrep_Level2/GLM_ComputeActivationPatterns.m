function [AnalysisParameters] = GLM_ComputeActivationPatterns(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
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
            [filePaths_glm,~,glmName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','ActivationPatterns','TitleTextName','Select Analysis Parameters:');
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
    RunContrast=AnalysisParameters.RunContrast;
    ConfoundNamesByRun=AnalysisParameters.ConfoundNamesByRun;
    ConfoundNamesBySubj=AnalysisParameters.ConfoundNamesBySubj;
    try
        fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    catch
        fmriprep_table_name=[];
    end
    EventNames=AnalysisParameters.EventNames;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    AM2=AnalysisParameters.AM2;
    ConditionNames=AnalysisParameters.ConditionNames;
    BaseTCName=AnalysisParameters.BaseTCName;  
    AnalysisName=AnalysisParameters.AnalysisName;
else
    AnalysisParameters=struct;
    ConditionNames=[];
    ConfoundNamesByRun=[];
    ConfoundNamesBySubj=[];
    EventNames=[];
    AM2Event_Names=[];
    TimeCourseNames=[];
    AnalysisParameters.AnalysisName=AnalysisName;
end
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
AnalysisParameters.fmriprep_table_name=fmriprep_table_name;
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

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
AnalysisParameters.UseAfni=UseAfni;
AnalysisParameters.Use3DReml=Use3DReml;
AnalysisParameters.AfniWorkDir=AfniWorkDir; 

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
ExperimentsDir=strrep(ExperimentsDir,'\','/');
%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run','Session','Group'},'Perform analysis by subject or by run:',SingleSelect);
end
bySes=0;
bySS=0;
byGroup=0;
byRun=0;
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=find(fmriprep_table.numRuns_bySub)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
elseif strcmpi(SubjectOrRun,'Group')
    byGroup=1;
    useIndicies=find(fmriprep_table.numRuns_byGroup)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_byGroup;    
elseif strcmpi(SubjectOrRun,'Session')
    bySes=1;
    useIndicies=find(fmriprep_table.numRuns_bySes)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;      
else    
    byRun=1;
    useIndicies=[1:TotalRuns];
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ActivationPatterns';
%Allows you to set name for this particular analysis
  
%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
[filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName',AnalysisName);
filePaths_beh_events=filePaths_beh_events.(beh_eventsName);

%% Select variables to include in GLM analysis
% Select confound regressors to include.
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
        if bySS == 1 
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
AnalysisParameters.ConfoundNamesBySubj=ConfoundNamesBySubj;
AnalysisParameters.ConfoundNamesByRun=ConfoundNamesByRun;
%Select Events and Timecourses
EventNamesAll=[];
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
AnalysisParameters.AM2=AM2;
AnalysisParameters.EventNames=EventNames;
AnalysisParameters.TimeCourseNames=TimeCourseNames;
AnalysisParameters.AM2Event_Names=AM2Event_Names;

%% Set contrast parameters if desired.
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
if isempty(AnalysisParameters.AnalysisName)
    tempName=strrep(['ActivationPattern_',afniSuffix],'__','_');
else
   tempName= AnalysisParameters.AnalysisName;
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

AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.BaseTCName=BaseTCName;
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
save('Temp_GLM_ComputeAP_Params', 'AnalysisParameters');
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
    
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1 && fmriprep_table.numRuns(dataInd,1)==fmriprep_table.numRuns_bySub(dataInd,1)
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end

    descript1='desc-activation_patterns'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    SavePrefix=[ExperimentsDir,SaveDir,'WholeBrain/']; 
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);   
    end
    
    SaveNames=[SavePrefix,SaveName,'.mat'];
    if exist(SaveNames,'file')~=0 && Overwrite==0
        disp(['Skipping-- all files exist: ',SaveNames]);
        continue       
    end
    %% Determine number of runs
    if bySS==1
        numRuns=fmriprep_table.numRuns_bySub(dataInd,1);
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
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
%         if ~isempty(fmriprep_table.Error{loadInd,1})
%             continue
%         end    
        %% set load paths and variable names
        %pull load paths
        LoadPath_Events=filePaths_beh_events{loadInd,1};
        LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
        LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
        %Pull base timecourse data. If it doesn't exist, skip run.
        try
            variableInfo = who('-file', LoadPath_BaseTC);
            if any(contains(variableInfo,'boldTCs'))
                TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                BaseTC{count,1}=TempLoadData.boldTCs;
                UseAfniConfounds=1;
            elseif any(contains(variableInfo,'Resids'))   
                TempLoadData = load(LoadPath_BaseTC,'Resids','brain_mask');
                BaseTC{count,1}=TempLoadData.Resids;
                UseAfniConfounds=0;
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
        
        count=count+1;
    end  
   
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    %% Run analysis here!!
    try
      [tVals,brain_mask,RegressorNames,~,DesignMatrix,AfniInfo,bVals] = GetActivationPatterns(...
            Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
            BaseTC,brain_mask,ConditionNames,...
            'UseAfni',UseAfni,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'Use3DReml',Use3DReml,...
            'UseAfniConfounds',UseAfniConfounds,...
            'TR',TR,...
            'TrialNum',TrialNum,...
            'parGLM',1,...
            'ContrastNames',[],...
            'Normalize','zscore',...
            'ConfoundTCsBySubj',Use_ConfoundTCsBySubj,...
            'ResampleSizes',RunDur,...
            'AM2',AM2,...
            'AM2_Events',Use_AM2Events,...
            'TimecourseShift',0,...
            'ParcellationMask',brain_mask{1,1},...
            'ResampleSizesBoldTC',[],...
            'NormWithinRun',[],...
            'NormAcrossRun','zscore',...
            'Contrasts',cMat,...
            'BatchSize',BatchSize,...
            'ResampleToFit','Y');   
    catch
        delete([AfniWorkDir,'*']);
        disp(['Error - ',SaveNames]);
        continue
    end
    save(SaveNames,'tVals','bVals','brain_mask','RegressorNames','DesignMatrix','AfniInfo','AnalysisParameters');     
    toc
end
end