function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeSubSamples(fmriprep_table,ExperimentsDir,varargin)
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
%Run analyses on TimeCourses input here (instead of selecting previously
%saved TimeCourses. Format is table (matched to fmriprep_table)
%of cells containing the timecourse to use. 
%[TimeCourses] = VariableSetter('TimeCourses',[],varargin);

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='subsample';
%Allows you to set name for this particular analysis

%% Compile filepaths for input files for the analysis
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
[filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
filePaths_beh_events=filePaths_beh_events.(beh_eventsName);

%% Select variables to include in GLM analysis
% Select confound regressors to include.
ConfoundNames=[];
EventNamesAll=[];
for aNum=1:TotalRuns
    if ~isempty(filePaths_ConfoundTCs{aNum,1})            
        load(filePaths_ConfoundTCs{aNum,1},'ConfoundTCs');            
        ConfoundNames=unique([ConfoundNames;ConfoundTCs.Properties.VariableNames(:)]);
    end    
end 
for aNum=1:TotalRuns
    if ~isempty(filePaths_beh_events{aNum,1})            
        load(filePaths_beh_events{aNum,1},'beh_events');            
        EventNamesAll=unique([EventNamesAll;beh_events.Properties.VariableNames(:)]);
    end    
end   
[SubsampleParams,SubSampleName] = uiTCSubSample([EventNamesAll(:);ConfoundNames(:)]);
TCName=SubsampleParams.SplitVar;
NumSubs=SubsampleParams.NumSubs;
SubSampleSize=SubsampleParams.SubSampleSize;

splitConditionNames=uiNameSelect([{'None'};EventNamesAll(:)],'Select conditions to subsample within',0);
if ~iscell(splitConditionNames)
    splitConditionNames={splitConditionNames};
end
if any(ismember(EventNamesAll,TCName))
    filePaths_TC=filePaths_beh_events;
    loadVarName='beh_events';
else
    filePaths_TC=filePaths_ConfoundTCs;
    loadVarName='ConfoundTCs';
end

if isempty(AnalysisName)
    AnalysisName=uiEnterName([AnalysisType,'_',SubSampleName],['Enter name for ',AnalysisType,newline,'analysis below:']);
end


%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

for rep=1:SubsampleParams.NumReps
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
    iniPercentComplete=0; %Used to display progress
    for dataInd=useIndicies
        tic
        %% Display progress
        PercentComplete=round((dataInd/TotalRuns)*100);
        if PercentComplete>iniPercentComplete
            disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName,': rep ',num2str(rep),' of ',num2str(SubsampleParams.NumReps)]);
            iniPercentComplete=PercentComplete;
        end
        %% Set save directory and save name
        SaveDir=[fmriprep_table.matDir{dataInd,1},AnalysisType,'/',AnalysisName,'/'];
        SaveName = fmriprep_table.preproc_bold{dataInd,1};
        SaveName=strrep(SaveName,'.nii','');
        SaveName=strrep(SaveName,'.gz','');
        if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
            SaveName=strrep(SaveName,'_run-01','');
        end
        descript1='desc-subsample'; %set file description name
        SaveName=strrep(SaveName,'desc-preproc_bold',descript1);

        %% If analysis does not involves parcellations:
        SavePrefix=[ExperimentsDir,SaveDir];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
        end    
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0 && rep==1 && SubsampleParams.SaveAppend==0
            disp(['Skipping-- file exists: ',SaveName]);
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
        UseSubSampleTC=cell(1);
        UseSubSampleConditions=cell(1);
        TrialNum=cell(1);
        RunDur=cell(1);
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
            LoadPath_TC=filePaths_TC{loadInd,1};
            LoadPath_events=filePaths_beh_events{loadInd,1};
            LoadPath_confoundTC=filePaths_ConfoundTCs{loadInd,1};
            %Pull base timecourse data. If it doesn't exist, skip run.        
            try
                TempLoadData = load(LoadPath_TC,loadVarName);
                UseSubSampleTC{count,1}=TempLoadData.(loadVarName).(TCName);
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_TC]);
                continue
            end 
            try
                TempLoadData = load(LoadPath_confoundTC,'ConfoundTCs');
                RunDur{count,1}=height(TempLoadData.ConfoundTCs);
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_confoundTC]);
                continue
            end         
            try
                if strcmpi(splitConditionNames{1,1},'None')
                    TempLoadData = load(LoadPath_events,'beh_events');
                    UseSubSampleConditions{count,1}=ones(length(UseSubSampleTC{count,1}),1);
                    TrialNum{count,1}=TempLoadData.beh_events.TrialNum;
                else
                    TempLoadData = load(LoadPath_events,'beh_events');
                    UseSubSampleConditions{count,1}=table2array(TempLoadData.beh_events(:,splitConditionNames));
                    TrialNum{count,1}=TempLoadData.beh_events.TrialNum;
                end
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_events]);
                continue
            end        
            count=count+1;
        end  

        if count==1
            disp(['Skipping subject-- no input files or variables exist']);
            continue
        end  
        %% Run analysis here!!
        SubConditions=cell(NumSubs,1);
        SubConditionNames=cell(NumSubs,1);
        SubAssigns=cell(NumSubs,1);
        SubData=cell(NumSubs,1);
        SubAssignsByTime=cell(NumSubs,1);
        skip=zeros(NumSubs,1);        
        parfor subNum=1:NumSubs
            try
                [tempSubConditions,tempSubConditionNames,tempSubAssigns,tempSubData,tempSubAssignsByTime] = ComputeSubsampleConditions(UseSubSampleConditions{1,1},splitConditionNames,UseSubSampleTC{1,1},subNum,NumSubs,'TrialNum',TrialNum{1,1},'SubSampleSize',SubSampleSize);
            catch
                skip(subNum,1)=1;
                continue                
            end
            SubConditions{subNum,1}=tempSubConditions;
            SubConditionNames{subNum,1}=tempSubConditionNames;
            SubAssigns{subNum,1}=tempSubAssigns;
            SubData{subNum,1}=[tempSubData.groupMeans,tempSubData.Error];
            SubAssignsByTime{subNum,1}=tempSubAssignsByTime;
        end
        if any(skip == 1)
            disp(['Error! Skipping: ',SaveNames{1,1}]);
            continue
        end
            
        SubConditionNames=SubConditionNames{1,1};
        SubAssignsByTrial=squeeze(cell2nDMAT(SubAssigns));
        SubAssignsByTime=squeeze(cell2nDMAT(SubAssignsByTime));
        SubAssignsByVol=imresize(SubAssignsByTime,[RunDur{1,1},NumSubs],'nearest');
        SubConditions=squeeze(cell2nDMAT(SubConditions));
        
        if rep == 1 && SubsampleParams.SaveAppend==0
            save(SaveNames{1,1},'SubAssignsByTrial','SubAssignsByTime','SubData','SubConditions','SubAssignsByVol','SubConditionNames','SubsampleParams','TrialNum','UseSubSampleTC','RunDur');
        else
            AppendSave(SaveNames{1,1},'SubData',SubData,'Append');
            AppendSave(SaveNames{1,1},'SubAssignsByTrial',SubAssignsByTrial,'Append');
            AppendSave(SaveNames{1,1},'SubAssignsByTime',SubAssignsByTime,'Append'); 
            AppendSave(SaveNames{1,1},'SubAssignsByVol',SubAssignsByVol,'Append'); 
            AppendSave(SaveNames{1,1},'SubConditions',SubConditions,'Append'); 
        end
     
        toc
    end
end 