function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeBehaviorTCs(fmriprep_table,ExperimentsDir,varargin)
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
[SubjectOrRun] = VariableSetter('SubjectOrRun',['Run'],varargin);

[IncludeResponse] = VariableSetter('IncludeResponse',[],varargin);
[InputDuration] = VariableSetter('InputDuration',[],varargin);
[mseq] = VariableSetter('mseq',[],varargin);

[UseParameters] = VariableSetter('UseParameters',[],varargin);

[PullParams] = VariableSetter('PullParams',[],varargin);
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='beh';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['events'],['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%% Compile filepaths for input files for the analysis
  
[filePaths_beh_raw,~,beh_raw_name] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','raw');
filePaths_beh_raw=filePaths_beh_raw.(beh_raw_name);

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

if isempty(IncludeResponse)
    SingleSelect=1; %Allows only a single value to be selected.
    [IncludeResponse] = uiNameSelect({'Yes','No'},'Include responses in behavior?',SingleSelect);
end
if strcmpi(IncludeResponse,'Yes')
    IncludeResponse=1;
else
    IncludeResponse=0;
end

if isempty(InputDuration)
    SingleSelect=1; %Allows only a single value to be selected.
    [InputDuration] = uiNameSelect({'Yes','No'},'Input duration?',SingleSelect);
end
if strcmpi(InputDuration,'Yes')
    InputDuration=1;
else
    InputDuration=0;
end
if isempty(mseq)
    SingleSelect=1; %Allows only a single value to be selected.
    [mseq] = uiNameSelect({'Yes','No'},'Include m-sequence timecourses?',SingleSelect);
end
if strcmpi(mseq,'Yes')
    mseq=2;
else
    mseq=0;
end

StimulusVariableName=[];
ResponseVariableName=[];
OnsetVariableName=[];
DurationVariableName=[];
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
redo=1;
while redo==1
    redoInd=[];
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
        descript1='desc-beh_events'; %set file description name
        SaveName=strrep(SaveName,'desc-preproc_bold',descript1);

        %% If analysis does not involves parcellations:
        SavePrefix=[ExperimentsDir,SaveDir,'/'];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
        end   
        SavePrefixClean=strrep(SavePrefix,['/',AnalysisName,'/'],'/clean/');    
        if ~exist(SavePrefixClean,'file')
            mkdir(SavePrefixClean);
        end        
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        SaveNames{2,1}=[SavePrefixClean,strrep(SaveName,descript1,'desc-beh_clean'),'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
            disp(['Skipping-- file exists: ',SaveName]);
            continue
        elseif exist(SaveNames{1,1},'file')~=0 && Overwrite==1 
            if isempty(PullParams)
                [PullParams] = uiNameSelect({'Yes','No'},'Pull exisiting parameters?',SingleSelect);
                if strcmpi(PullParams,'Yes')
                    PullParams=1;
                else
                    PullParams=0;
                end
            end
            if PullParams==1    
                load(SaveNames{1,1},'UseParameters')
            end
        end

        %% Initialize input data for loading
        if bySS==1
            numRuns=fmriprep_table.numRuns(dataInd,1);
        else
            numRuns=1;
        end
        beh_raw=[];
        %% load input data 
        % If by Subject, iterate through runs and place data in cell

        %skip previous errors
        if ~isempty(fmriprep_table.Error{dataInd,1})
            continue
        end    
        %% set load paths and variable names
        %pull load paths
        LoadPath_beh_raw=filePaths_beh_raw{dataInd,1};

        %Pull base timecourse data. If it doesn't exist, skip run.
        try
            TempLoadData = load(LoadPath_beh_raw,'beh_raw');
            beh_raw=TempLoadData.beh_raw;
        catch
            disp(['Skipping run-- input file or variable doesnt exist']);
            continue
        end 

        %% Run analysis here!!
        if any(ismember(fmriprep_table.Properties.VariableNames,'fmriDur'))
            fmriDur=fmriprep_table.fmriDur{dataInd,1};
        else
            fmriDur=[];
        end
        if isempty(StimulusVariableName)
            StimulusVariableName=uiNameSelect(beh_raw.Properties.VariableNames,'Select stimulus variable header.',0);
        end    
        StimVars=[];
        if iscell(StimulusVariableName)
            for i= 1:size(StimulusVariableName,1)
                StimVars=[StimVars,beh_raw.(StimulusVariableName{i,1})];
            end
            StimVars=join(StimVars,'_',2);
        else
            StimVars=beh_raw.(StimulusVariableName);
        end
        if isempty(OnsetVariableName)
            OnsetVariableName=uiNameSelect(beh_raw.Properties.VariableNames,'Select onset variable header.',1);
        end  
        Onsets=beh_raw.(OnsetVariableName);
        TrialDuration=[];
        if isempty(DurationVariableName) && InputDuration==1
            DurationVariableName=uiNameSelect(beh_raw.Properties.VariableNames,'Select duration variable header.',1);
        end   
        if InputDuration==1
            TrialDuration=beh_raw.(DurationVariableName);
        end
        
        Responses=[];
        if isempty(ResponseVariableName) && IncludeResponse==1
            ResponseVariableName=uiNameSelect(beh_raw.Properties.VariableNames,'Select response variable header.',1);
        end    
        if IncludeResponse==1
            Responses=beh_raw.(ResponseVariableName);
        end
        
        try
            [ beh_clean,beh_events,UseParameters ] = ExtractEvents_General(StimVars,Onsets,...
                'Responses',Responses,...
                'TrialDuration',TrialDuration,...
                'UseParameters',UseParameters,...
                'mseq',mseq,...
                'fmriDur',fmriDur);
        catch
            redoInd=[redoInd,dataInd];
            continue
        end
            save(SaveNames{1,1},'beh_events','UseParameters'); 
            save(SaveNames{2,1},'beh_clean','UseParameters'); 

        toc
    end
    if isempty(redoInd)
        redo=0;
    else
        useIndicies=redoInd;
        UseParameters=[];
        StimulusVariableName=[];
        ResponseVariableName=[];
        OnsetVariableName=[];
    end
end
end