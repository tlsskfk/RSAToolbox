function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = GLM_ComputeSubSampleRSMs(fmriprep_table,ExperimentsDir,varargin)
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
[AM2] = VariableSetter('AM2',0,varargin);
%Run GLM Contrast (default is [] triggering a UI selection)
[RunContrast] = VariableSetter('RunContrast',0,varargin);
%Set batch size for glm. set for optimal tradeoff between speed and memory.
[BatchSize] = VariableSetter('BatchSize',1000,varargin);
% Set Parcellation to run analysis on
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
% Set type of similarity measure for RSMs
[SimType] = VariableSetter('SimType',[],varargin);
% Compute Confound RSMs (Cell containing names of confound RSMs to compute)
[ComputeConfoundRSMs] = VariableSetter('ComputeConfoundRSMs',[],varargin);
[UseStat] = VariableSetter('UseStat',[],varargin); %T of B
[SubSamplesName] = VariableSetter('SubSamplesName',[],varargin); 
[BaseTCName] = VariableSetter('SubSamplesName',[],varargin); 
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','subsample_RSMs','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_RSMs.Properties.VariableNames;
        catch
            filePaths_RSMs=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_RSMs)
            for j = 1:height(filePaths_RSMs)
                try
                    load(filePaths_RSMs.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
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
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    SimType=AnalysisParameters.SimType;
    UseStat=AnalysisParameters.UseStat;
    ComputeConfoundRSMs=AnalysisParameters.ComputeConfoundRSMs;
    gmMask=AnalysisParameters.gmMask;
    ConfoundNamesByRun=AnalysisParameters.ConfoundNamesByRun;
    ConfoundNamesBySubj=AnalysisParameters.ConfoundNamesBySubj;
    EventNames=AnalysisParameters.EventNames;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    AM2=AnalysisParameters.AM2;
    ParcelNames=AnalysisParameters.ParcelNames;
    SubSamplesName=AnalysisParameters.SubSampleName;
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
AnalysisParameters.SubjectOrRun=SubjectOrRun;
if isempty(SimType)
    SimType = uiNameSelect({'corrcoef','meanSim','euclidean','squaredeuclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'},'Select similarity/distance measure to use.',1);     
end
AnalysisParameters.SimType=SimType;
if isempty(UseStat)
    UseStat = uiNameSelect({'T','B'},'Activation pattern value to use.',1);     
end
AnalysisParameters.UseStat=UseStat;
%% Compile filepaths for input files for the analysis
[filePaths_SubSamples,~,SubSamplesName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','subsample','AnalysisName',SubSamplesName,'TitleTextName','Select subsample type for RSM:');
filePaths_SubSamples=filePaths_SubSamples.(SubSamplesName);
AnalysisParameters.SubSampleName=SubSamplesName;
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='subsample_RSMs';
AnalysisParameters.AnalysisType=AnalysisType;
%% Specify which Confound RSMs (if any) to compute.
Compute_CondReg_CorrMat=0;
Compute_WM_RSM=0;
Compute_CSF_RSM=0;
Compute_GM_RSM=0;
Compute_MeanAct_RSM=0;
ConfoundAffix=['_cl'];
if isempty(ComputeConfoundRSMs)
    ComputeConfoundRSMs = uiNameSelect({'None','CSF','WM','GM','CondReg_CorrMat','MeanAct'},'Select Confound RSMs to compute.');     
end
AnalysisParameters.ComputeConfoundRSMs=ComputeConfoundRSMs;
if any(ismember(ComputeConfoundRSMs,'CondReg_CorrMat'))
    Compute_CondReg_CorrMat=1;
    ConfoundAffix=[ConfoundAffix,'Rg'];
end
if any(ismember(ComputeConfoundRSMs,'WM'))
    Compute_WM_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Wm'];
end
if any(ismember(ComputeConfoundRSMs,'CSF'))
    Compute_CSF_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Cf'];
end
if any(ismember(ComputeConfoundRSMs,'GM'))
    Compute_GM_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Gm'];
end
if any(ismember(ComputeConfoundRSMs,'MeanAct'))
    Compute_MeanAct_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Mn'];
end
if strcmpi(ConfoundAffix,'_cl')
    ConfoundAffix=['_'];
end

%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
AnalysisParameters.BaseTCName=BaseTCName;
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
[filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
filePaths_beh_events=filePaths_beh_events.(beh_eventsName);
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
        gmAffix=['_gm',num2str4filename(gmMask,2)];
    end
end
if strcmpi(gmMask,'NA_EMPTY')
    gmMask=[];
end
AnalysisParameters.gmMask=gmMask;
if ~isempty(gmMask) || Compute_WM_RSM || Compute_CSF_RSM || Compute_GM_RSM
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end
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
    if AM2==1  && ~isempty(TimeCourseNames) && ~strcmpi(AM2Event_Names,'NA_EMPTY')
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
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    tempName=[AnalysisType,gmAffix,ConfoundAffix,'_',SubSamplesName];
    tempName=strrep(tempName,'subsample_RSMs','subRSMs');
    tempName=strrep(tempName,'SubSample','sub');
    tempName=strrep(tempName,'subsample','sub');
    tempName=strrep(tempName,'_sub_','_');
    tempName=strrep(tempName,'_cl_','_');
    tempName=strrep(tempName,'__','_');
    AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
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
    ParcelNames=uiNameSelect([{'Searchlight'};ParcelNames],'Select parcellations to run.');
elseif iscell(ParcelNames)
    ParcelNames=ParcelNames(:);
end
AnalysisParameters.ParcelNames=ParcelNames;
numParcels=length(ParcelNames);
SearchlightParams=[];
if contains(ParcelNames,'Searchlight')
    SearchlightParams=uiSearchlight;
end

%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];
AllParcels=cell(numParcels,2);
if ~any(ismember(ParcelNames,'Searchlight')) && Compute_WM_RSM==0 && Compute_CSF_RSM==0 && Compute_GM_RSM==0
    for i = 1:numParcels
        load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
        AllParcels{i,1}=UseMask;
        AllParcels{i,2}=UseLabels;
        if i ==1            
            analysisMask=single(UseMask>0);
        else
            analysisMask=analysisMask+single(UseMask>0);
        end
        analysisMask=single(analysisMask>0);
    end
end
BaseAnalysisMask=analysisMask;
for i = 1:numParcels
    if ~strcmpi(ParcelNames{i,1},'Searchlight')
        load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
        AllParcels{i,1}=UseMask;
        AllParcels{i,2}=UseLabels;
    else
        AllParcels{i,1}=[];
        AllParcels{i,2}=[]; 
    end
end

Parcels=cell(numParcels,1);
for i = 1:numParcels
    if ~strcmpi(ParcelNames{i,1},'Searchlight')
        Parcels{i,1}=load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
    end
end

if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end
AnalysisParameters.SearchlightParams=SearchlightParams;

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    analysisMask=BaseAnalysisMask;
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
    descript1='desc-subrsms'; %set file description name
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
    Use_ConfoundEvents=cell(1);
    Use_TimeCourses=cell(1);
    Use_AM2Events=cell(1);
    Use_gmMask=[];    
    Use_SubSamples=cell(1);
    SubConditionNames=cell(1);
    AllSubData=cell(1);
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
        LoadPath_subsamples=filePaths_SubSamples{loadInd,1};
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
                    ConfoundEvents=table2array(Events);
                    ConfoundEventNames=Events.Properties.VariableNames(:);
                else
                    ConfoundEvents=[];
                    ConfoundEventNames=[];
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
        end
        %Pull subsample events
        
        if ~isempty(LoadPath_subsamples)
            try
                TempLoadData = load(LoadPath_subsamples,'SubConditions','SubConditionNames','SubData');               
                numReps=size(TempLoadData.SubConditions,4);
                numSubs=size(TempLoadData.SubConditions,3);                
                for r = 1:numReps
                    for s = 1:numSubs
                        Use_SubSamples{r,s}{count,1}=[TempLoadData.SubConditions(:,:,s,r),ConfoundEvents];
                    end
                end
                if count==1
                    ConditionNames=TempLoadData.SubConditionNames(1:end-1,1);
                    SubConditionNames=[TempLoadData.SubConditionNames;ConfoundEventNames]; 
                end
                SubData=cell2mat(TempLoadData.SubData);
                AllSubData{1,count}=SubData(:,1);
            catch
                disp(['No events-- input file or variable doesnt exist: ',LoadPath_subsamples]);
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
        
        if Compute_WM_RSM == 1
            try
                TempLoadData = load(filePaths_gmMask{loadInd,1},'WM_probseg');
                Use_wmMask=single(TempLoadData.WM_probseg > 0.99);
            catch
                disp(['White matter mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                Use_wmMask=[];
            end
        else
            Use_wmMask=[];
        end 
        
        if Compute_CSF_RSM == 1
            try
                TempLoadData = load(filePaths_gmMask{loadInd,1},'CSF_probseg');
                Use_csfMask=single(TempLoadData.CSF_probseg > 0.9);
            catch
                disp(['CSF mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                Use_csfMask=[];
            end
        else
            Use_csfMask=[];
        end     
         if Compute_GM_RSM == 1
            try
                TempLoadData = load(filePaths_gmMask{loadInd,1},'GM_probseg');
                Use_invgmMask=single(TempLoadData.GM_probseg < 0.01);
            catch
                disp(['invGM mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                Use_invgmMask=[];
            end
        else
            Use_invgmMask=[];
        end        
        count=count+1;
    end  
   
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    if Compute_CSF_RSM == 0 && Compute_WM_RSM == 0 && Compute_GM_RSM == 0
        if isempty(analysisMask)
            analysisMask=Use_gmMask;
        else
            analysisMask=analysisMask.*Use_gmMask;
        end
    else
        analysisMask=(brain_mask{1,1}*0)+1;        
    end
    %% Run analysis here!!
    
    for repNum=1:numReps
        AllRSMs=cell(numParcels,numSubs);
        all_rsm_mask=cell(numParcels,numSubs);
        temp_Use_Subsamples=Use_SubSamples(repNum,:);
        AllSubData=cell2mat(AllSubData);
        subsample_data=nanmean(AllSubData,2);         
        try
            parfor subNum=1:numSubs                  
                Use_Events=[temp_Use_Subsamples{1,subNum},repmat({SubConditionNames},[size(temp_Use_Subsamples{1,subNum},1),1])];
                [AllRSMs(:,subNum),all_rsm_mask(:,subNum)] = GetParcellationRSMs(...
                    Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
                    BaseTC,brain_mask,ConditionNames,Parcels,...
                    'csfMask',Use_csfMask,...
                    'wmMask',Use_wmMask,...
                    'invgmMask',Use_invgmMask,...
                    'Compute_MeanAct_RSM',Compute_MeanAct_RSM,...
                    'gmMask',Use_gmMask,...
                    'SimType',SimType,...
                    'OutputConfoundRSMs',0,...
                    'parGLM',0,...
                    'Compute_CondReg_CorrMat',Compute_CondReg_CorrMat,...
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
                    'ResampleToFit','Y',...
                    'UseStat',UseStat);                             
            end
        catch
            disp(['Error! Skipping: ', SaveNames{1,1}])
            continue            
        end
        AllRSMs=AllRSMs';
        all_rsm_mask=all_rsm_mask';
        for i = 1:numParcels
            RSMs=single(cell2nDMAT(AllRSMs(:,i)));
            rsm_mask=all_rsm_mask{1,i};
            if repNum == 1
                save(SaveNames{i,1},'RSMs','rsm_mask','AnalysisParameters','subsample_data');
            else
                AppendSave(SaveNames{i,1},'RSMs',RSMs,'Append');
                AppendSave(SaveNames{i,1},'subsample_data',subsample_data,'Append');
            end
        end            
    end
    toc
end
end 

function residRSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs)
    numConfoundRSMs=size(ConfoundRSMs,1);
    RSMSize=size(RSMs,1);
    RSMs = mat2uppertriuvectormat(RSMs);
    numRSMs=size(RSMs,2);
    indMat=ones(numRSMs,numConfoundRSMs);
    residRSMs=RSMs*0;
    RegConstant=ones(size(RSMs,1),1);
    for i = 1:numConfoundRSMs
        ConfoundRSMs{i,1}=mat2uppertriuvectormat(ConfoundRSMs{i,1});
        if size(ConfoundRSMs{i,1},2)==numRSMs
            indMat(:,i)=[1:numRSMs]';
        end
    end
    parfor i = 1:numRSMs
        tempRSMs=RSMs(:,i);
        tempConfoundRSMs=RegConstant;
        for j = 1:numConfoundRSMs
            tempConfoundRSMs=[tempConfoundRSMs,ConfoundRSMs{j,1}(:,indMat(i,j))];
        end
        [residRSMs(:,i)] = FastOLSRegress_Resids(tempRSMs,tempConfoundRSMs);
    end
    %[ CleanRSMs ] = vertRSM2SymRSM( residRSMs,RSMSize );    
end

