function [fmriprep_table] = FunctionTemplate(fmriprep_table,ExperimentsDir,varargin)
%Template function for data processing from the BIDsTable

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);

%For other variable inputs XX use VariableSetter function where 'Variable' is the
%text string indicating the variable name and DefaultVal is the value
%assigned to Variable if nothing is specified.
[Variable] = VariableSetter('Variable',DefaultVal,varargin);


%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='AnalysisType';

%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(genDateString,['Enter name for ',AnalysisType,newline,'analysis below:']);
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

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%Create Table column to record the location of the output of this
%analysis. May become obsolete.
if sum(ismember(fmriprep_table.Properties.VariableNames,[AnalysisType,'_',AnalysisName]))==0
    fmriprep_table.([AnalysisType,'_',AnalysisName])=cell(TotalRuns,1);
end

%% Search a mat file for names to select from
ConfoundNamesAll=[];
for aNum=1:numAnalyses
    load([ExperimentsDir,fmriprep_table.LocalPath_mat{aNum,1},fmriprep_table.ConfoundTCs{aNum,1}],'regressor_names');
    ConfoundNamesAll=unique([ConfoundNamesAll;regressor_names(:)]);
end
%Alternatively use this function:
[outNames] = BIDsDirSearch(fmriprep_table,ExperimentsDir,varargin);

%Select a subset of items from a larger list
[EventNames] = uiNameSelect([UniqueEventNames],'Select events or blocks to include:');


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
    ParcelNames=uiNameSelect([ParcelNames;{'Searchlight'}],'Select parcellations to run.');
else
    ParcelNames=Parcellation(:);
end
numParcels=length(ParcelNames);

%If running searchlight analysis, set searchlight parameters
if contains(ParcelNames,'Searchlight')
    SearchlightParams=uiSearchlight;
end

%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
ParcelSum_Mask=[];
if sum(ismember(ParcelNames,'Searchlight'))==0
    for i = 1:numParcels
        load(['Parcellations/',ParcelNames{i,1}],'UseMask');
        if i ==1            
            ParcelSum_Mask=single(UseMask>0);
        else
            ParcelSum_Mask=ParcelSum_Mask+single(UseMask>0);
        end
        ParcelSum_Mask=single(ParcelSum_Mask>0);
    end
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
    if fmriprep_table.numRuns{dataInd,1}>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-XXXXX'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    fmriprep_table.([AnalysisType,'_',AnalysisName]){dataInd,1}={SaveDir,SaveName};
    
    %% If analysis involves parcellations:
    SavePrefix=[ExperimentsDir,SaveDir,ParcelNames{parcelNum,1},'/'];    
    SaveNames=cell(numParcels,1);
    RunParcel=ones(numParcels,1);     
    for parcelNum=1:numParcels
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
        numRuns=fmriprep_table.numRuns{dataInd,1};
    else
        numRuns=1;
    end
    count=1;
    InputData=cell(1);
    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
        if ~isempty(fmriprep_table.Error{loadInd,1})
            continue
        end    
        %% set load paths and variable names
        LoadVars={'VarName1','VarName2'};
        LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.LocalPath_mat{loadInd,1},fmriprep_table.(LoadVar1){dataInd,1}];
        LoadPaths{1,2}=[ExperimentsDir,fmriprep_table.LocalPath_mat{loadInd,1},fmriprep_table.(LoadVar2){dataInd,1}];
        for loadNum=1:length(LoadPaths)
            try
                TempLoadData = load(LoadPaths{1,loadNum},LoadVars{1,loadNum});
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,loadNum},newline,LoadPaths{1,loadNum}]);
                continue
            end       
            InputData{count,loadNum}=TempLoadData.LoadVar; 
            InputData{count,loadNum}=TempLoadData.LoadVar;
        end
        count=count+1;
    end  
    
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ',LoadVars{1,1},newline,LoadPaths{1,1}]);
        continue
    end  
    
    %% Run analysis here!!
    [Condition_tVals,brain_mask,AllRegressorNames] = GetActivationPatterns(...
        AllEvents,AllTimeCourses,AllConfoundTCs,...
        AllboldTCs,AllboldMasks,ConditionNames,...
        'parGLM',1,...
        'ContrastNames',cNames,...
        'Normalize','zscore',...
        'ConfoundTCsBySubj',AllConfoundTCsBySubj,...
        'ResampleSizes',ResampDur,...
        'AM2',AM2,...
        'AM2_Events',All_AM2_Events,...
        'TimecourseShift',TimecourseShift,...
        'ParcellationMask',ParcelSum_Mask,...
        'ResampleSizesBoldTC',[],...
        'NormWithinRun',[],...
        'NormAcrossRun','zscore',...
        'Contrasts',cMat,...
        'BatchSize',BatchSize,...
        'ResampleToFit','Y'); 
    
    for parcelNum=1:numParcels
        if SaveBySS{parcelNum,2}==0
            continue
        end
        load(['Parcellations/',ParcelNames{parcelNum,1}],'UseMask');
        ParcellationVector=UseMask(brain_mask~=0); 
        rsm_mask=single(brain_mask~=0).*UseMask;
        [ RSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
        save(SaveBySS{parcelNum,1},'RSMs','rsm_mask','ConditionNames','AllRegressorNames');   
    end    
    toc
end
end 