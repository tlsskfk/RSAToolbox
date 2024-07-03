function [AnalysisParameters] = GLM_ComputeRSMs(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
end
%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Run AM2 regression for timecourses (default is 1 for yes)
[LoadActivationPatterns] = VariableSetter('LoadActivationPatterns',[],varargin);
[ActivationPatternName] = VariableSetter('ActivationPatternName',[],varargin);
[AM2] = VariableSetter('AM2',1,varargin);
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
% Voxel Normalization
[VoxNorm] = VariableSetter('VoxNorm',[],varargin);
% Compute Confound RSMs (Cell containing names of confound RSMs to compute)
[ComputeConfoundRSMs] = VariableSetter('ComputeConfoundRSMs',[],varargin);
SearchlightName = VariableSetter('SearchlightName',[],varargin);
SearchlightParams = VariableSetter('SearchlightParams',[],varargin);
[BaseTCName] = VariableSetter('SubSamplesName',[],varargin);
[UseStat] = VariableSetter('UseStat',[],varargin); %T of B
[UseAfni] = VariableSetter('UseAfni',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
Use3DReml = VariableSetter('Use3DReml',0,varargin);
if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','RSMs','TitleTextName','Select Analysis Parameters:');
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
    try
        VoxNorm=AnalysisParameters.VoxNorm;  
        UseAfni=AnalysisParameters.UseAfni;
        AfniWorkDir=AnalysisParameters.AfniWorkDir;
        Use3DReml=AnalysisParameters.Use3DReml;
    catch
        VoxNorm=[];
        UseAfni=0;
        AfniWorkDir=[];
        Use3DReml=0;
    end
    SearchlightName=AnalysisParameters.SearchlightName;
    SearchlightParams=AnalysisParameters.SearchlightParams;
    LoadActivationPatterns=AnalysisParameters.LoadActivationPatterns;
    ComputeConfoundRSMs=AnalysisParameters.ComputeConfoundRSMs;
    gmMask=AnalysisParameters.gmMask;
    ActivationPatternName=AnalysisParameters.ActivationPatternName;
    ConfoundNamesByRun=AnalysisParameters.ConfoundNamesByRun;
    ConfoundNamesBySubj=AnalysisParameters.ConfoundNamesBySubj;
    EventNames=AnalysisParameters.EventNames;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    AM2=AnalysisParameters.AM2;
    ConditionNames=AnalysisParameters.ConditionNames;
    ParcelNames=AnalysisParameters.ParcelNames;
    BaseTCName=AnalysisParameters.BaseTCName;   
    ParcelType=AnalysisParameters.ParcelType;  
    
else
    AnalysisParameters=struct;
    ConditionNames=[];
    ConfoundNamesByRun=[];
    ConfoundNamesBySubj=[];
    EventNames=[];
    AM2Event_Names=[];
    TimeCourseNames=[];
    cNames=[];
end

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
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
AnalysisType='RSMs';
AnalysisParameters.AnalysisType=AnalysisType;
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
%% Specify which Confound RSMs (if any) to compute.
Compute_CondReg_CorrMat=0;
Compute_WM_RSM=0;
Compute_CSF_RSM=0;
Compute_GM_RSM=0;
Compute_MeanAct_RSM=0;
ConfoundAffix=['_cl'];

if isempty(SimType)
    SimType = uiNameSelect({'corrcoef','meanSim','euclidean','squaredeuclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'},'Select similarity/distance measure to use.',1);     
end
AnalysisParameters.SimType=SimType;
if isempty(VoxNorm)
    VoxNorm = uiNameSelect({'None','mean','Z'},'Select voxel normalization type.',1);     
end
AnalysisParameters.VoxNorm=VoxNorm;
if strcmpi(VoxNorm,'None')
    VoxNormAffix='';
else
    VoxNormAffix=['_Vox',VoxNorm];
end

if isempty(UseStat)
    UseStat = uiNameSelect({'T','B'},'Activation pattern value to use.',1);     
end
AnalysisParameters.UseStat=UseStat;
if isempty(ComputeConfoundRSMs)
    ComputeConfoundRSMs = uiNameSelect({'None','CSF','WM','GM','CondReg_CorrMat','MeanAct'},'Select Confound RSMs to compute.');     
end
if strcmpi(ComputeConfoundRSMs,'NA_EMPTY')
    ComputeConfoundRSMs=[];
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
if any(ismember(ComputeConfoundRSMs,'GM'))
    Compute_GM_RSM=1;
    ConfoundAffix=[ConfoundAffix,'Gm'];
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
if LoadActivationPatterns==0
    [filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
    filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
    AnalysisParameters.BaseTCName=BaseTCName;
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
AnalysisParameters.ConditionNames=ConditionNames;
if LoadActivationPatterns==1
    [filePaths_ActivationPattern,~,ActivationPatternName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ActivationPatterns','AnalysisName',ActivationPatternName,'TitleTextName','Select activation patterns for RSA:');
    filePaths_ActivationPattern=filePaths_ActivationPattern.WholeBrain;
end
AnalysisParameters.ActivationPatternName=ActivationPatternName;
%Allows you to set name for this particular analysis
if LoadActivationPatterns==1
    UseActivationPatternName=strrep(ActivationPatternName,'ActivationPattern_','RSMs_');
    UseActivationPatternName=strrep(UseActivationPatternName,'ActivationPatterns_','RSMs_');
    tempName=[UseActivationPatternName,gmAffix,ConfoundAffix,'_',SimType,UseStat,VoxNormAffix];
else
    tempName=[AnalysisType,gmAffix,ConfoundAffix,'_',SimType,UseStat,afniSuffix,VoxNormAffix];
end
if isempty(AnalysisName)
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
        ParcelNames=uiNameSelect([{'Searchlight'};ParcelNames],'Select parcellations to run.');
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
        [tempParcelFileNames] = sub2runFilePath(fmriprep_table.numRuns_bySub,tempParcelFileNames);
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
LoadSLName=0;
RunSearchlight=0;
if contains(ParcelNames,'Searchlight')
    RunSearchlight=1;
    if isempty(SearchlightName) && isempty(SearchlightParams)
        SingleSelect=1; %Allows only a single value to be selected.
        [LoadSLName] = uiNameSelect({'Yes','No'},'Load searchlight info?:',SingleSelect);
        if strcmpi(LoadSLName,'Yes')
            LoadSLName=1;
            SearchlightParams=1;
            [filePaths_SearchlightInfo,~,SearchlightName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SearchlightInfo','AnalysisName',SearchlightName,'TitleTextName','Select searchlight info:');
            filePaths_SearchlightInfo=filePaths_SearchlightInfo.(SearchlightName);
        else
            filePaths_SearchlightInfo=[];
            LoadSLName=0;         
            SingleSelect=1; %Allows only a single value to be selected.
            [slShape] = uiNameSelect({'Sphere','Cube'},'Select searchlight shape:',SingleSelect);
            slRadiusName=uiEnterName('3','Enter searchlight radius');
            slRadius=str2num(slRadiusName);   
            if strcmpi(slShape,'Sphere')
                maxVoxels=MaxSphereSLSize(slRadius,2);
            else
                maxVoxels=(slRadius*2+1)^3;
            end
            slVoxelThreshold=uiEnterName(num2str(round(maxVoxels*0.75)),['Enter searchlight radius ',newline,'Max # voxels: ',num2str(maxVoxels)]);
            SearchlightName=uiEnterName(['SearchlightInfo_Shape-',slShape,'_Rad-',slRadiusName,'_Thresh-',slVoxelThreshold],['Enter searchlight name:']);
            slVoxelThreshold=str2num(slVoxelThreshold);
            SearchlightParams.slShape=slShape;
            SearchlightParams.slRadius=slRadius;
            SearchlightParams.slVoxelThreshold=slVoxelThreshold;
        end
        
    end
else
    SearchlightName='';
end
SearchlightName=strrep(SearchlightName,'SearchlightInfo_','Searchlight_');
AnalysisParameters.SearchlightName=SearchlightName;
AnalysisParameters.SearchlightParams=SearchlightParams;
Parcels=cell(numParcels,1);
if strcmpi(ParcelType,'GroupParcel')
    for i = 1:numParcels
        if ~strcmpi(ParcelNames{i,1},'Searchlight')
            Parcels{i,1}=load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
        else
            Parcels{i,1}='Searchlight';
            ParcelNames{i,1}=SearchlightName;
        end
    end

    %Create mask that excludes voxels that aren't included in any parcellation.
    %Helpful to save time and minimize computation.
    analysisMask=[];

    if ~any(ismember(ParcelNames,SearchlightName)) && Compute_WM_RSM==0 && Compute_CSF_RSM==0 && Compute_GM_RSM==0
        for i = 1:numParcels
            UseMask=Parcels{i,1}.UseMask;
            if i ==1            
                analysisMask=single(UseMask>0);
            else
                try
                    analysisMask=analysisMask+single(UseMask>0);
                catch
                    disp('oops')
                    continue
                end
            end
            analysisMask=single(analysisMask>0);
        end
    end
else
   analysisMask=[];
end
BaseAnalysisMask=analysisMask;
if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end
MeanRSMs=cell(numParcels,2);
MeanZRSMs=cell(numParcels,2);
AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.SearchlightParams=SearchlightParams;
AnalysisParameters.BaseTCName=BaseTCName;
save('Temp_GLM_ComputeRSMs_Params', 'AnalysisParameters');
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    try
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
            runNum=fmriprep_table.run(dataInd,1);
            SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
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
        Use_ActivationPatterns=[];
        Use_Events=cell(1);
        Use_TimeCourses=cell(1);
        Use_AM2Events=cell(1);
        Use_SearchlightInfo=[];
        Use_gmMask=[];    
        TR=fmriprep_table.TR(dataInd,1); 
        TrialNum=cell(1);     
        TrialOnsetTimes=cell(1);    
        %% load input data 
        % If by Subject, iterate through runs and place data in cell
        if LoadActivationPatterns==1
            LoadPath_ActivationPatterns=filePaths_ActivationPattern{dataInd,1};
            if strcmpi(UseStat,'T')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'tVals','brain_mask','AnalysisParameters');
                    Use_ActivationPatterns.ActVals=Use_ActivationPatterns.tVals;
                    Use_ActivationPatterns.tVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end
            elseif strcmpi(UseStat,'B')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'bVals','brain_mask','AnalysisParameters');
                    Use_ActivationPatterns.ActVals=Use_ActivationPatterns.bVals;
                    Use_ActivationPatterns.bVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end 
            end
            brain_mask=repmat(brain_mask,[numRuns,1]);
            AnalysisParameters.ActivationPatterns=Use_ActivationPatterns.AnalysisParameters;
        end   
        if LoadSLName==1
            LoadPath_SearchlightInfo=filePaths_SearchlightInfo{dataInd,1};
            Use_SearchlightInfo = load(LoadPath_SearchlightInfo,'SLcoords','SLinds','SLnumVoxels','AnalysisParameters');
            Use_SearchlightInfo.slShape=Use_SearchlightInfo.AnalysisParameters.slShape;
            Use_SearchlightInfo.slRadius=Use_SearchlightInfo.AnalysisParameters.slRadius;
            Use_SearchlightInfo.slVoxelThreshold=Use_SearchlightInfo.AnalysisParameters.slVoxelThreshold;
        else
            if RunSearchlight == 1
                Use_SearchlightInfo=SearchlightParams;
                Use_SearchlightInfo.SLcoords=[];
                Use_SearchlightInfo.SLinds=[];
                Use_SearchlightInfo.SLnumVoxels=[];
            end
        end
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
                try
                    LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
                
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
                    continue
                    disp(['invGM mask error-- mask not applied']);
                    continue
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
        if ~strcmpi(ParcelType,'GroupParcel') 
            tempParcelPaths=table2cell(AllParcelTable(dataInd,:));
            Parcels=cell(size(tempParcelPaths,2),1);
            for j = 1:size(tempParcelPaths,2)
                Parcels{j,1}=load(tempParcelPaths{1,j},'UseMask','UseLabels');          
                if unique(Parcels{j,1}.UseMask(:)) == 0
                    disp(['Skipping ',tempParcelPaths{1,j},' -- ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                    continue
                end
            end
        end
        
        %% Run analysis here!!
        [All_RSMs,All_rsm_masks,All_ConfoundRSMs,AnalysisParameters.RegressorNames,AnalysisParameters.AfniInfo,SearchlightResults] = GetParcellationRSMs(...
            Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
            BaseTC,brain_mask,ConditionNames,Parcels,...
            'SearchlightInfo',Use_SearchlightInfo,...
            'ActivationPatterns',Use_ActivationPatterns,...
            'UseAfni',UseAfni,...
            'Use3DReml',Use3DReml,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'TR',TR,...
            'TrialNum',TrialNum,...            
            'csfMask',Use_csfMask,...
            'wmMask',Use_wmMask,...
            'invgmMask',Use_invgmMask,...
            'Compute_MeanAct_RSM',Compute_MeanAct_RSM,...
            'gmMask',Use_gmMask,...
            'SimType',SimType,...
            'OutputConfoundRSMs',0,...
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
            'ResampleToFit','Y',...
            'VoxNorm',VoxNorm,...
            'UseStat',UseStat);   
        for i = 1:numParcels
            if strcmpi(ParcelNames{i,1},'Searchlight')
                RSMs=single(All_RSMs{i,1});
                rsm_mask=All_rsm_masks{i,1};                
                save(SaveNames{i,1},'RSMs','rsm_mask','SearchlightResults','AnalysisParameters');
            else 
                RSMs=All_RSMs{i,1};
                rsm_mask=All_rsm_masks{i,1};
                if ~isempty(All_ConfoundRSMs)
                    confoundRSMs=All_ConfoundRSMs{i,1};
                    save(SaveNames{i,1},'RSMs','rsm_mask','confoundRSMs','AnalysisParameters');   
                else
                    save(SaveNames{i,1},'RSMs','rsm_mask','AnalysisParameters');   
                end  
            end
        end 
        toc
    catch
       disp('Error!');
    end
end 
for parcelNum = 1:numParcels
    parcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];  
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
    end
    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);     
    end
    if contains(parcelName,'Searchlight')
        [all_masks] = CompileND(ExperimentsDir,fmriprep_table,...
            'LoadVarName','rsm_mask',...
            'LoadVarFormat','Array',...
            'DataType','ByParcellation',...
            'AnalysisType',AnalysisType,...
            'AnalysisName',AnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'ParcelName',parcelName);
        OverlapMap=sum(all_masks,4);
        filterVec=OverlapMap(:)==0;
        OverlapVec=OverlapMap(:);
        OverlapVec(filterVec,:)=[];
        slMask=OverlapMap>0;
        numSS=size(all_masks,4);
        slPaths=BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',AnalysisName);
        slPaths=slPaths.(parcelName);
        All_RSMs=[];
        count=1;
        for s=1:size(slPaths,1)
            try
                slData=load(slPaths{s,1},'rsm_mask','RSMs');
            catch
                continue
            end
            if isempty(All_RSMs)
                All_RSMs=nan(size(slData.RSMs,1),length(OverlapVec),numSS,'single');
            end
            tempInd=slData.rsm_mask(:);
            tempInd(filterVec,:)=[];
            All_RSMs(:,tempInd>0,count)=slData.RSMs;
            count=count+1;
        end  
        slRSMs=nanmean(All_RSMs,3);
        for s=1:size(All_RSMs,3)
            All_RSMs(:,:,s)=nan_zscore(All_RSMs(:,:,s));
        end
        slZRSMs=nanmean(All_RSMs,3);
        save([GroupAnalysisDir,'slRSMs.mat'],'slRSMs','slZRSMs','OverlapVec','slMask');
        save([GroupAnalysisDir,'AnalysisParameters.mat'],'AnalysisParameters');        
    else    
        [All_RSMs] = CompileND(ExperimentsDir,fmriprep_table,...
            'LoadVarName','RSMs',...
            'LoadVarFormat','Array',...
            'DataType','ByParcellation',...
            'AnalysisType',AnalysisType,...
            'AnalysisName',AnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'ParcelName',parcelName); 
        if ndims(All_RSMs) == 2 && length(Parcels{parcelNum,1}.UseLabels) == 1
            All_RSMs=reshape(All_RSMs,[size(All_RSMs,1),1,size(All_RSMs,2)]);
        end       
        All_RSMs_Z=All_RSMs*0;

        for n= 1:size(All_RSMs,3)
            All_RSMs_Z(:,:,n)=nan_zscore(All_RSMs(:,:,n));
        end
        MeanRSM=nanmean(All_RSMs,3);
        ParcelRSMs=array2table(MeanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanZRSM=nanmean(All_RSMs_Z,3);
        ParcelZRSMs=array2table(MeanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanCleanRSM=MeanRSM*0;
        MeanCleanZRSM=MeanCleanRSM;
        [~,RF]=FastRCA(All_RSMs);
        RF=RF';
        numROIs=size(All_RSMs,2);
        for roi=1:numROIs
            MeanCleanRSM(:,roi)=nanmean(squeeze(All_RSMs(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
            MeanCleanZRSM(:,roi)=nanmean(squeeze(All_RSMs_Z(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
        end   
        ParcelRSMs_p50=array2table(MeanCleanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');
        ParcelZRSMs_p50=array2table(MeanCleanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');
        save([GroupAnalysisDir,'ParcelRSMs.mat'],'ParcelRSMs','ParcelZRSMs','ParcelRSMs_p50','ParcelZRSMs_p50');
        save([GroupAnalysisDir,'AnalysisParameters.mat'],'AnalysisParameters');       
    end
end
end