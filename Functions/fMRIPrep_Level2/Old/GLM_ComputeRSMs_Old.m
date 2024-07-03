function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = GLM_ComputeRSMs(fmriprep_table,ExperimentsDir,varargin)
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
[RunContrast] = VariableSetter('RunContrast',0,varargin);
%Set batch size for glm. set for optimal tradeoff between speed and memory.
[BatchSize] = VariableSetter('BatchSize',1000,varargin);
% Set Parcellation to run analysis on
[Parcellation] = VariableSetter('Parcellation',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
% Set type of similarity measure for RSMs
[SimType] = VariableSetter('SimType','corrcoef',varargin);
% Compute Confound RSMs (Cell containing names of confound RSMs to compute)
[ComputeConfoundRSMs] = VariableSetter('ComputeConfoundRSMs',[],varargin);


AnalysisParams=struct;
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
AnalysisType='RSMs';

%% Specify which Confound RSMs (if any) to compute.
Compute_CondReg_CorrMat=0;
Compute_WM_RSM=0;
Compute_CSF_RSM=0;
Compute_MeanAct_RSM=0;
ConfoundAffix=['_cl'];
if isempty(ComputeConfoundRSMs)
    ComputeConfoundRSMs = uiNameSelect({'None','CSF','WM','CondReg_CorrMat','MeanAct'},'Select Confound RSMs to compute.');     
end

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
if any(ismember(ComputeConfoundRSMs,'MeanAct'))
    Compute_MeanAct_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Mn'];
end
if strcmpi(ConfoundAffix,'cl')
    ConfoundAffix=['_clNone'];
end

%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
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

if ~isempty(gmMask) || Compute_WM_RSM || Compute_CSF_RSM
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end

%% Select variables to include in GLM analysis
% Select confound regressors to include.
ConfoundNames=[];
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

%Select Events and Timecourses
EventNamesAll=[];
EventNames=[];
TimeCourseNames=[];
AM2Event_Names=[];
if any(cellfun(@isempty,filePaths_beh_events)==0)
    for aNum=1:size(filePaths_beh_events,1)
        if ~isempty(filePaths_beh_events{aNum,1})            
            load(filePaths_beh_events{aNum,1},'beh_events');            
            EventNamesAll=unique([EventNamesAll;beh_events.Properties.VariableNames(:)]);
        end    
    end
    EventNames=uiNameSelect(EventNamesAll,'Select events to include');
    TimeCourseNames=uiNameSelect(EventNamesAll,'Select timecourses to include');
    if AM2==1  && ~isempty(TimeCourseNames)
        for TCNum=1:length(TimeCourseNames)
            [AM2Event_Name] = uiNameSelect([EventNamesAll],['Select AM2 Event for ',TimeCourseNames{TCNum,1},':']);
            AM2Event_Names=[AM2Event_Names;AM2Event_Name];
        end    
    end    
else
    AM2=0;
end

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
ConditionNames = uiNameSelect([EventNames(:);TimeCourseNames(:);cNames(:);unique(AM2Event_Names(:));ConfoundNamesByRun(:);ConfoundNamesBySubj(:)],'Select conditions/contrasts to save:');

%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName([AnalysisType,gmAffix,ConfoundAffix,'_',genDateString],['Enter name for ',AnalysisType,newline,'analysis below:']);
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
if ischar(Parcellation)
    ParcelNames={Parcellation};
elseif ~iscell(Parcellation)
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
else
    ParcelNames=Parcellation(:);
end
numParcels=length(ParcelNames);
SearchlightParams=[];
if contains(ParcelNames,'Searchlight')
    SearchlightParams=uiSearchlight;
end
%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];
if ~any(ismember(ParcelNames,'Searchlight')) && Compute_WM_RSM==0 && Compute_CSF_RSM==0
    for i = 1:numParcels
        load(['Parcellations/',ParcelNames{i,1}],'UseMask');
        if i ==1            
            analysisMask=single(UseMask>0);
        else
            analysisMask=analysisMask+single(UseMask>0);
        end
        analysisMask=single(analysisMask>0);
    end
end

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end

AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.Parcellation=ParcelNames;
AnalysisParams.gmMask=gmMask;
AnalysisParams.ConfoundRSMs=ComputeConfoundRSMs;
AnalysisParams.SearchlightParams=SearchlightParams;
AnalysisParams.ConditionNames=ConditionNames;
AnalysisParams.BaseTCName=BaseTCName;
AnalysisParams.SimType=SimType;

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    try
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
        end
        descript1='desc-rsms'; %set file description name
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
        Use_TimeCourses=cell(1);
        Use_AM2Events=cell(1);
        Use_gmMask=[];    
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
                    Use_wmMask=single(TempLoadData.WM_probseg > 0.9);
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
            count=count+1;
        end  

        if count==1
            disp(['Skipping subject-- no input files or variables exist']);
            continue
        end  
        if Compute_CSF_RSM == 0 && Compute_WM_RSM == 0
            if isempty(analysisMask)
                analysisMask=Use_gmMask;
            else
                analysisMask=analysisMask.*Use_gmMask;
            end
        else
            analysisMask=(brain_mask{1,1}*0)+1;        
        end
        %% Run analysis here!!
        [Condition_tVals,brain_mask,AnalysisParams.RegressorNames,CondReg_CorrMat] = GetActivationPatterns(...
            Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
            BaseTC,brain_mask,ConditionNames,...
            'parGLM',1,...
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
            'ResampleToFit','Y');   
        brainSize=size(brain_mask);
        ConfoundRSMs=cell(1);
        ConfoundRSMLabels=cell(1);
        if Compute_CondReg_CorrMat==1
            ConfoundRSMs=cat(1,ConfoundRSMs,{CondReg_CorrMat});
            ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'CondReg_CorrMat'});
        end    
        if Compute_CSF_RSM == 1 && ~isempty(Use_csfMask)
            ParcellationVector=Use_csfMask(brain_mask~=0); 
            [ CSF_RSM ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
            ConfoundRSMs=cat(1,ConfoundRSMs,{CSF_RSM});
            ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'CSF_RSM'});       
        else
            CSF_RSM=[];
        end
        if Compute_WM_RSM == 1 && ~isempty(Use_wmMask)
            ParcellationVector=Use_wmMask(brain_mask~=0); 
            [ WM_RSM ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
            ConfoundRSMs=cat(1,ConfoundRSMs,{WM_RSM});
            ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'WM_RSM'});        
        else
            WM_RSM=[];
        end    
        if ~isempty(ConfoundRSMs)
            if isempty(ConfoundRSMs{1,1})
                ConfoundRSMs(1,:)=[];
                ConfoundRSMLabels(1,:)=[];
            end
            numConfoundRSMs=size(ConfoundRSMs,1);
        end

        for i = 1:numParcels
            if strcmpi(ParcelNames{i,1},'Searchlight')
                brainVals=Condition_tVals;
                save(SaveNames{i,1},'brainVals','brain_mask','AnalysisParams');
            else           
                load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
                if ~ismember(brainSize,size(UseMask),'rows')
                    UseMask=imresize3(UseMask,brainSize,'nearest');
                end    
                UseLabels=UseLabels(:);
                numLabels=length(UseLabels);
                gmVector=Use_gmMask(brain_mask~=0);
                ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
                rsm_mask=single(brain_mask~=0).*UseMask.*Use_gmMask;
                [ RSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
                if Compute_MeanAct_RSM==1
                    [ MeanActRSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,'meanSim');
                    if isempty(ConfoundRSMs)
                        ConfoundRSMs{1,1}=MeanActRSMs;
                    else
                        ConfoundRSMs{numConfoundRSMs+1,1}=MeanActRSMs;
                        ConfoundRSMLabels{numConfoundRSMs+1,1}='MeanActRSMs';
                    end
                end
                if ~isempty(ConfoundRSMs)
                    RSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs);
                end
                save(SaveNames{i,1},'RSMs','rsm_mask','AnalysisParams');   
            end
        end 
        toc
    catch
        disp('Error!');
    end
end 
end
function CleanRSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs)
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
    [ CleanRSMs ] = vertRSM2SymRSM( residRSMs,RSMSize );    
end

