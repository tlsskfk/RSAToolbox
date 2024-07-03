function [AnalysisParameters] = Compute2SplitVarParcellationRCA(fmriprep_table,ExperimentsDir,varargin)
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
[LoadActivationPatterns] = VariableSetter('LoadActivationPatterns',[],varargin);
[ActivationPatternName] = VariableSetter('ActivationPatternName',[],varargin);
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);
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
%Run GLM Contrast (default is [] triggering a UI selection)
[SplitAcrossRun] = VariableSetter('SplitAcrossRun',[],varargin);
[DownSample] = VariableSetter('DownSample',[],varargin);
%Subject or run level analysis. Will prompt request.
[DownSampleN] = VariableSetter('DownSampleN',[],varargin);
%Subject or run level analysis. Will prompt request.
% Voxel Normalization
[VoxNorm] = VariableSetter('VoxNorm',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
[DownSampleReps] = VariableSetter('DownSampleReps',100,varargin);
%Subject or run level analysis. Will prompt request.
[RunPermute] = VariableSetter('RunPermute',[],varargin);
%Subject or run level analysis. Will prompt request.
[PermuteReps] = VariableSetter('PermuteReps',[],varargin);
%Subject or run level analysis. Will prompt request.
[OutputPermDists] = VariableSetter('OutputPermDists',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[SkipRCA] = VariableSetter('SkipRCA',0,varargin);
zNormRSM = VariableSetter( 'zNormRSM',[],varargin);
zNormRCA = VariableSetter( 'zNormRCA',[],varargin);
SearchlightName = VariableSetter('SearchlightName',[],varargin);
SearchlightParams = VariableSetter('SearchlightParams',[],varargin);
%Input predefined analysis parameters analysis parameters
[UseAfni] = VariableSetter('UseAfni',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
Use3DReml = VariableSetter('Use3DReml',1,varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_SplitRCA,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','SplitVarRSMs','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_SplitRCA.Properties.VariableNames;
        catch
            filePaths_SplitRCA=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_SplitRCA)
            for j = 1:height(filePaths_SplitRCA)
                try
                    load(filePaths_SplitRCA.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
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
    LoadActivationPatterns=AnalysisParameters.LoadActivationPatterns;
    ActivationPatternName=AnalysisParameters.ActivationPatternName;  
    FigureInfo=AnalysisParameters.FigureInfo;
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    SimType=AnalysisParameters.SimType;
    UseStat=AnalysisParameters.UseStat;
    DownSample=AnalysisParameters.DownSample;
    DownSampleN=AnalysisParameters.DownSampleN;
    DownSampleReps=AnalysisParameters.DownSampleReps;
    RunPermute=AnalysisParameters.RunPermute;
    PermuteReps=AnalysisParameters.PermuteReps;
    OutputPermDists=AnalysisParameters.OutputPermDists;    
    ComputeConfoundRSMs=AnalysisParameters.ComputeConfoundRSMs;
    zNormRSM=AnalysisParameters.zNormRSM;
    zNormRCA=AnalysisParameters.zNormRCA;
    SearchlightName=AnalysisParameters.SearchlightName;
    SearchlightParams=AnalysisParameters.SearchlightParams;
    gmMask=AnalysisParameters.gmMask;
    UseAfni=AnalysisParameters.UseAfni;
    AfniWorkDir=AnalysisParameters.AfniWorkDir;   
    Use3DReml=AnalysisParameters.Use3DReml;
    ConfoundNamesByRun=AnalysisParameters.ConfoundNamesByRun;
    ConfoundNamesBySubj=AnalysisParameters.ConfoundNamesBySubj;
    EventNames=AnalysisParameters.EventNames;
    TimeCourseNames=AnalysisParameters.TimeCourseNames;
    AM2Event_Names=AnalysisParameters.AM2Event_Names;
    VoxNorm=AnalysisParameters.VoxNorm;
    AM2=AnalysisParameters.AM2;
    ConditionNames=AnalysisParameters.ConditionNames;
    ParcelNames=AnalysisParameters.ParcelNames;
    SplitAcrossRun=AnalysisParameters.SplitAcrossRun;
    %NumReps=AnalysisParameters.NumReps;
    ConfoundEventNames=AnalysisParameters.ConfoundEventNames;
    SplitVar=AnalysisParameters.SplitVar;
    BaseTCName=AnalysisParameters.BaseTCName;
    try
        fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    catch
        fmriprep_table_name=[];
    end
else
    AnalysisParameters=struct;
    ConditionNames=[];
    ConfoundNamesByRun=[];
    ConfoundNamesBySubj=[];
    EventNames=[];
    AM2Event_Names=[];
    TimeCourseNames=[];
    ConfoundEventNames=[];
    SplitVar=[];
    BaseTCName=[];
    cNames=[];    
end
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants
if isempty(fmriprep_table_name)
[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

% Set figure output info
if isempty(FigureInfo)
    SingleSelect=0; %Allows only a single value to be selected.
    [FigureInfo] = uiNameSelect({'None','MakeFigures','ComputeReliability','MakeBrainMaps'},'Select summary figures to make: ',SingleSelect);    
end
AnalysisParameters.FigureInfo=FigureInfo;
MakeFigs=0;
ComputeReliability=0;
MakeBrainMaps=0;
if any(contains(FigureInfo,'MakeFigures'))
    MakeFigs=1;
end
if any(contains(FigureInfo,'ComputeReliability'))
    ComputeReliability=1;
end
if any(contains(FigureInfo,'MakeBrainMaps'))
    MakeBrainMaps=1;
end

if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
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
    SplitAcrossRun=0;
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;
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
if LoadActivationPatterns==1
    [filePaths_ActivationPattern,~,ActivationPatternName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ActivationPattern2SplitVar','AnalysisName',ActivationPatternName,'TitleTextName','Select Activation Patterns for RSA:');
    filePaths_ActivationPattern=filePaths_ActivationPattern.WholeBrain;
end
AnalysisParameters.ActivationPatternName=ActivationPatternName;
if isempty(zNormRSM)
    SingleSelect=1; %Allows only a single value to be selected.
    [zNormRSM] = uiNameSelect({'Yes','No'},'Z normalize RSMs for RCA?',SingleSelect);
    if strcmpi(zNormRSM,'Yes')
        zNormRSM=1;
    else
        zNormRSM=0;
    end    
end
if zNormRSM == 1
    zNormSuffix='_zRSM';
else
    zNormSuffix='';
end
AnalysisParameters.zNormRSM=zNormRSM;
if isempty(zNormRCA)
    SingleSelect=1; %Allows only a single value to be selected.
    [zNormRCA] = uiNameSelect({'Yes','No'},'Z normalize RCA values?',SingleSelect);
    if strcmpi(zNormRCA,'Yes')
        zNormRCA=1;
    else
        zNormRCA=0;
    end    
end
if zNormRCA == 1
    zNormRCASuffix='_zRCA';
else
    zNormRCASuffix='';
end
if isempty(VoxNorm)
    VoxNorm = uiNameSelect({'None','mean','Z'},'Select voxel normalization type.',1);     
end
AnalysisParameters.VoxNorm=VoxNorm;
if strcmpi(VoxNorm,'None')
    VoxNormAffix='';
else
    VoxNormAffix=['_Vox',VoxNorm];
end
AnalysisParameters.zNormRSM=zNormRSM;
AnalysisParameters.zNormRCA=zNormRCA;
if isempty(DownSample)
    SingleSelect=1; %Allows only a single value to be selected.
    [DownSample] = uiNameSelect({'Yes','No'},'Down sample RCA to match a group N:',SingleSelect);
    if strcmpi(DownSample,'Yes')
        DownSample=1;
    else
        DownSample=0;
    end
end
if DownSample==1 
    if isempty(DownSampleN)
        DownSampleName=uiEnterName('','Enter downsample N');
        DownSampleN=str2num(DownSampleName);
    else
        DownSampleName=num2str(DownSampleN);
    end
    if isempty(DownSampleReps)
        DownSampleReps=uiEnterName('100','Enter # downsample reps');
        DownSampleReps=str2num(DownSampleReps);
    end
        
    DownSampleSuffix=['_N',DownSampleName];
else
    DownSampleSuffix='';
end
AnalysisParameters.DownSample=DownSample;
AnalysisParameters.DownSampleN=DownSampleN;
AnalysisParameters.DownSampleReps=DownSampleReps;
if isempty(RunPermute)
    SingleSelect=1; %Allows only a single value to be selected.
    [RunPermute] = uiNameSelect({'Yes','No'},'Run permutation analysis:',SingleSelect);
    if strcmpi(RunPermute,'Yes') 
        RunPermute=1;
    else
        RunPermute=0;
    end    
end
if RunPermute==1 
    if isempty(PermuteReps)
        PermuteRepsName=uiEnterName('','Enter # permutation reps:');
        PermuteReps=str2num(PermuteRepsName);
    else
        PermuteRepsName=num2str(PermuteReps);
    end
    PermuteSuffix=['_perm',PermuteRepsName];
    if isempty(OutputPermDists)
        SingleSelect=1; %Allows only a single value to be selected.
        [OutputPermDists] = uiNameSelect({'Yes','No'},'Save permutations by ss/run:',SingleSelect);
        if strcmpi(OutputPermDists,'Yes')
            OutputPermDists=1;
        else
            OutputPermDists=0;
        end
    end    
else
    PermuteSuffix='';
end
AnalysisParameters.RunPermute=RunPermute;
AnalysisParameters.PermuteReps=PermuteReps;
AnalysisParameters.OutputPermDists=OutputPermDists;

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
afniSuffix='';
if LoadActivationPatterns==0
    if bySS==1
        if isempty(SplitAcrossRun)
            SingleSelect=1; %Allows only a single value to be selected.
            [SplitAcrossRun] = uiNameSelect({'Yes','No'},'Split across runs?',SingleSelect);
            if strcmpi(SplitAcrossRun,'Yes')
                SplitAcrossRun=1;
            else
                SplitAcrossRun=0;
            end
        end
    end

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
    AnalysisParameters.SplitAcrossRun=SplitAcrossRun;
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SplitVarRSMs';
AnalysisParameters.AnalysisType=AnalysisType;
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
if isempty(UseStat)
    UseStat = uiNameSelect({'T','B'},'Activation pattern value to use.',1);     
end
AnalysisParameters.UseStat=UseStat;
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
if strcmpi(ConfoundAffix,'cl')
    ConfoundAffix=['_clNone'];
end
    
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
        gmAffix=['_gm',num2str4filename(gmMask,2)];
    end
elseif isnumeric(gmMask)    
    gmAffix=['_gm',num2str4filename(gmMask,2)];
elseif strcmpi(gmMask,'NA_EMPTY')
    gmMask=[];
end
AnalysisParameters.gmMask=gmMask;
if ~isempty(gmMask) || Compute_WM_RSM || Compute_CSF_RSM || Compute_GM_RSM
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
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
        if isempty(EventNames) || isempty(TimeCourseNames) || isempty(AM2Event_Names) || isempty(SplitVar)
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
    if isempty(SplitVar)
        SplitVar=uiNameSelect([EventNamesAll(:);ConfoundNames(:)],'Select Split Var',1);
    end
    if iscell(SplitVar)
        SplitVar=SplitVar{1,1};
    end
    AnalysisParameters.AM2=AM2;
    AnalysisParameters.EventNames=EventNames;
    AnalysisParameters.TimeCourseNames=TimeCourseNames;
    AnalysisParameters.AM2Event_Names=AM2Event_Names;
    AnalysisParameters.SplitVar=SplitVar;
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

    if length(ConditionNames)<length(EventNames)
        ConfoundEventNames=EventNames(~ismember(EventNames,ConditionNames),:);
    end
    AnalysisParameters.ConfoundEventNames=ConfoundEventNames;
else
    BaseTCName=[];
    AM2=[];
    EventNames=[];
    TimeCourseNames=[];
    AM2Event_Names=[];
    SplitVar=[];
    ConfoundEventNames=[];
    ConditionNames=[];
    ConfoundNamesBySubj=[];
    ConfoundNamesByRun=[];
    AnalysisParameters.AM2=AM2;
    AnalysisParameters.EventNames=EventNames;
    AnalysisParameters.TimeCourseNames=TimeCourseNames;
    AnalysisParameters.AM2Event_Names=AM2Event_Names;
    AnalysisParameters.SplitVar=SplitVar;    
    AnalysisParameters.ConfoundEventNames=ConfoundEventNames;
    AnalysisParameters.ConditionNames=ConditionNames;
    AnalysisParameters.ConfoundNamesBySubj=ConfoundNamesBySubj;
    AnalysisParameters.ConfoundNamesByRun=ConfoundNamesByRun;    
end
%Allows you to set name for this particular analysis
tempName=[strrep(ActivationPatternName,'ActivationPattern_','RSMs_'),gmAffix,ConfoundAffix,'_',SimType,UseStat,VoxNormAffix];
if isempty(AnalysisName)
    tempName=strrep(tempName,'_cl_','_');
    tempName=strrep(tempName,'__','_');
    if DefaultName~=1
        AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
    else
        AnalysisName=tempName;
    end
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
        [tempParcelFileNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','CoordParcels','AnalysisName',ParcelNames{i,1});
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

AnalysisParameters.AnalysisName=AnalysisName;

AnalysisParameters.BaseTCName=BaseTCName;

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    %try
        tic
        analysisMask=BaseAnalysisMask;
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
            runNum=fmriprep_table.run(dataInd,1);
            SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
        end
        descript1='desc-SplitVarRSMs'; %set file description name
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
        Use_ConfoundEvents=cell(1);
        Use_TimeCourses=cell(1);
        Use_AM2Events=cell(1);
        TrialNum=cell(1);
        TrialOnsetTimes=cell(1); 
        Use_ActivationPatterns=[];
        Use_gmMask=[];  
        UseSplitVarTC=cell(1);
        TR=fmriprep_table.TR(dataInd,1);        
        %% load input data 
        % If by Subject, iterate through runs and place data in cell
        if LoadActivationPatterns==1
            LoadPath_ActivationPatterns=filePaths_ActivationPattern{dataInd,1};
            if strcmpi(UseStat,'T')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'Split1_tVals','Split2_tVals','brain_mask','DesignMatrix','AfniInfo','SplitData','AnalysisParameters');
                    Use_ActivationPatterns.Split1_ActVals=Use_ActivationPatterns.Split1_tVals;
                    Use_ActivationPatterns.Split2_ActVals=Use_ActivationPatterns.Split2_tVals;
                    Use_ActivationPatterns.Split1_tVals=[];
                    Use_ActivationPatterns.Split2_tVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end
            elseif strcmpi(UseStat,'B')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'Split1_bVals','Split2_bVals','brain_mask','DesignMatrix','AfniInfo','SplitData','AnalysisParameters');
                    Use_ActivationPatterns.Split1_ActVals=Use_ActivationPatterns.Split1_bVals;
                    Use_ActivationPatterns.Split2_ActVals=Use_ActivationPatterns.Split2_bVals;
                    Use_ActivationPatterns.Split1_bVals=[];
                    Use_ActivationPatterns.Split2_bVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end 
            end
            brain_mask=repmat(brain_mask,[numRuns,1]);
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
            else
                Use_SearchlightInfo=[];
            end
        end        
        for run=1:numRuns
            loadInd=dataInd+run-1;
            %skip previous errors
            if ~isempty(fmriprep_table.Error{loadInd,1})
                continue
            end    
            if LoadActivationPatterns==0 
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
                try
                    TempLoadData = load(LoadPath_Events,'beh_events');     
                    TrialNum{count,1}=TempLoadData.beh_events.TrialNum; 
                    TrialOnsetTimes{count,1}=TempLoadData.beh_events.TrialOnsetTime; 
                catch
                   disp(['Event load error-- skipping',LoadPath_Events]); 
                   continue
                end

                try
                    TempLoadData = load(LoadPath_Events,'beh_events');     
                    UseSplitVarTC{count,1}=TempLoadData.beh_events.(SplitVar); 
                catch
                   disp(['Event load error-- skipping',LoadPath_Events,SplitVar]); 
                   continue
                end            
                if ~isempty(LoadPath_Events)
                    try
                        TempLoadData = load(LoadPath_Events,'beh_events');
                        if ~isempty(ConditionNames)
                            Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,ConditionNames)));
                            Use_Events{count,1}=table2array(Events);
                            Use_Events{count,2}=Events.Properties.VariableNames(:);
                        end 
                        if ~isempty(ConfoundEventNames)
                            Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,ConfoundEventNames)));
                            Use_ConfoundEvents{count,1}=table2array(Events);
                            Use_ConfoundEvents{count,2}=Events.Properties.VariableNames(:);
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
        if ~strcmpi(ParcelType,'GroupParcel') 
            tempParcelPaths=table2cell(AllParcelTable(dataInd,:));
            Parcels=cell(size(tempParcelPaths,2),1);
            for j = 1:size(tempParcelPaths,2)
                Parcels{j,1}=load(tempParcelPaths{1,j},'UseMask','UseLabels');          
                if unique(Parcels{j,1}.UseMask(:)) == 0
                    disp(['Skipping ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                    continue
                end
            end
        end        
        %% Run analysis here!!
        try 
        [All_RSMs,brain_mask,SplitData] = GetParcellation2SplitVarRCA(...
            Use_Events,Use_ConfoundEvents,Use_TimeCourses,Use_ConfoundTCsByRun,...
            BaseTC,brain_mask,ConditionNames,Parcels,TrialNum,UseSplitVarTC,...
            'ActivationPatterns',Use_ActivationPatterns,...  
            'SearchlightInfo',Use_SearchlightInfo,...
            'UseAfni',UseAfni,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'TR',TR,...
            'Use3DReml',Use3DReml,...
            'SplitAcrossRun',SplitAcrossRun,...
            'csfMask',Use_csfMask,...
            'wmMask',Use_wmMask,...
            'invgmMask',Use_invgmMask,...
            'Compute_MeanAct_RSM',Compute_MeanAct_RSM,...
            'gmMask',Use_gmMask,...
            'SimType',SimType,...
            'OutputConfoundRSMs',0,...
            'parGLM',1,...
            'Compute_CondReg_CorrMat',Compute_CondReg_CorrMat,...
            'ContrastNames',[],...
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
            'Contrasts',[],...
            'BatchSize',BatchSize,...
            'ResampleToFit','Y',...
            'UseStat',UseStat); 
        catch
            disp(['Error: ',SaveName])
            continue
        end
        rsm_mask=brain_mask;
        for i = 1:numParcels
            Split1_RSMs=single(All_RSMs{i,1});
            Split2_RSMs=single(All_RSMs{i,2});
            save(SaveNames{i,1},'Split1_RSMs','Split2_RSMs','SplitData','brain_mask','rsm_mask','AnalysisParameters');                  
        end 
        toc
%     catch
%         disp(['Unknown Error-- ',SaveNames{i,1}]);
%     end
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
            'LoadVarName','brain_mask',...
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
        slPaths=BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SplitVarRSMs','AnalysisName',AnalysisName);
        slPaths=slPaths.(parcelName);
        All_RSMs=[];
        count=1;
        for s=1:size(slPaths,1)
            try
                slData=load(slPaths{s,1},'brain_mask','Split1_RSMs');
            catch
                continue
            end
            if isempty(All_RSMs)
                All_RSMs=nan(size(slData.Split1_RSMs,1),length(OverlapVec),numSS,'single');
            end
            tempInd=slData.brain_mask(:);
            tempInd(filterVec,:)=[];
            All_RSMs(:,tempInd>0,count)=slData.Split1_RSMs;
            count=count+1;
        end  
        Split1_slRSMs=nanmean(All_RSMs,3);
        for s=1:size(All_RSMs,3)
            All_RSMs(:,:,s)=nan_zscore(All_RSMs(:,:,s));
        end
        Split1_slZRSMs=nanmean(All_RSMs,3);
        
        All_RSMs=[];
        count=1;
        for s=1:size(slPaths,1)
            try
                slData=load(slPaths{s,1},'brain_mask','Split2_RSMs');
            catch
                continue
            end
            if isempty(All_RSMs)
                All_RSMs=nan(size(slData.Split2_RSMs,1),length(OverlapVec),numSS,'single');
            end
            tempInd=slData.brain_mask(:);
            tempInd(filterVec,:)=[];
            All_RSMs(:,tempInd>0,count)=slData.Split2_RSMs;
            count=count+1;
        end  
        Split2_slRSMs=nanmean(All_RSMs,3);
        for s=1:size(All_RSMs,3)
            All_RSMs(:,:,s)=nan_zscore(All_RSMs(:,:,s));
        end
        Split2_slZRSMs=nanmean(All_RSMs,3);
        
        save([GroupAnalysisDir,'slRSMs.mat'],'Split1_slRSMs','Split1_slZRSMs','Split2_slRSMs','Split2_slZRSMs','OverlapVec','slMask');
        save([GroupAnalysisDir,'AnalysisParameters.mat'],'AnalysisParameters');        
    else    
        [All_RSMs] = CompileND(ExperimentsDir,fmriprep_table,...
            'LoadVarName','Split1_RSMs',...
            'LoadVarFormat','Array',...
            'DataType','ByParcellation',...
            'AnalysisType',AnalysisType,...
            'AnalysisName',AnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'ParcelName',parcelName); 
        All_RSMs_Z=All_RSMs*0;
        for n= 1:size(All_RSMs,3)
            All_RSMs_Z(:,:,n)=nan_zscore(All_RSMs(:,:,n));
        end
        MeanRSM=nanmean(All_RSMs,3);
        Split1_ParcelRSMs=array2table(MeanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanZRSM=nanmean(All_RSMs_Z,3);
        Split1_ParcelZRSMs=array2table(MeanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanCleanRSM=MeanRSM*0;
        MeanCleanZRSM=MeanCleanRSM;
        [~,RF]=FastRCA(All_RSMs);
        RF=RF';
        numROIs=size(All_RSMs,2);
        for roi=1:numROIs
            MeanCleanRSM(:,roi)=nanmean(squeeze(All_RSMs(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
            MeanCleanZRSM(:,roi)=nanmean(squeeze(All_RSMs_Z(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
        end   
        Split1_ParcelRSMs_p50=array2table(MeanCleanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');
        Split1_ParcelZRSMs_p50=array2table(MeanCleanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');
        
        [All_RSMs] = CompileND(ExperimentsDir,fmriprep_table,...
            'LoadVarName','Split2_RSMs',...
            'LoadVarFormat','Array',...
            'DataType','ByParcellation',...
            'AnalysisType',AnalysisType,...
            'AnalysisName',AnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'ParcelName',parcelName); 
        All_RSMs_Z=All_RSMs*0;
        for n= 1:size(All_RSMs,3)
            All_RSMs_Z(:,:,n)=nan_zscore(All_RSMs(:,:,n));
        end
        MeanRSM=nanmean(All_RSMs,3);
        Split2_ParcelRSMs=array2table(MeanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanZRSM=nanmean(All_RSMs_Z,3);
        Split2_ParcelZRSMs=array2table(MeanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');

        MeanCleanRSM=MeanRSM*0;
        MeanCleanZRSM=MeanCleanRSM;
        [~,RF]=FastRCA(All_RSMs);
        RF=RF';
        numROIs=size(All_RSMs,2);
        for roi=1:numROIs
            MeanCleanRSM(:,roi)=nanmean(squeeze(All_RSMs(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
            MeanCleanZRSM(:,roi)=nanmean(squeeze(All_RSMs_Z(:,roi,RF(:,roi)>nanmedian(RF(:,roi)))),2);
        end   
        Split2_ParcelRSMs_p50=array2table(MeanCleanRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');
        Split2_ParcelZRSMs_p50=array2table(MeanCleanZRSM,'VariableNames',Parcels{parcelNum,1}.UseLabels');        
        
        
        save([GroupAnalysisDir,'ParcelRSMs.mat'],'Split1_ParcelRSMs','Split1_ParcelZRSMs','Split1_ParcelRSMs_p50','Split1_ParcelZRSMs_p50','Split2_ParcelRSMs','Split2_ParcelZRSMs','Split2_ParcelRSMs_p50','Split2_ParcelZRSMs_p50');
        save([GroupAnalysisDir,'AnalysisParameters.mat'],'AnalysisParameters');       
    end
end
if RunSearchlight == 1
    SkipRCA=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% Run RCA Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SkipRCA==0
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SplitVarRCA';
%Allows you to set name for this particular analysis

%% Compile filepaths for input files for the analysis
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SplitVarRSMs','AnalysisName',AnalysisName);
AnalysisName=strrep(RSMsName,'SplitVarRSMs_',['SplitVarRCA_',fmriprep_table_name,'_']);
AnalysisName=[AnalysisName,zNormSuffix,zNormRCASuffix,DownSampleSuffix,PermuteSuffix];

%% Identify parcellations availible and select ones to run analysis on.
% ParcelNames=filePaths_RSMs.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
%ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_RSMs=filePaths_RSMs(:,ParcelNames);
numParcels=length(ParcelNames);

%% define save filepaths from the load filepaths 
[LoadDirs_RSMs,LoadNames_RSMs] = filepath2DirAndFileName(filePaths_RSMs);
oldDesc='desc-';
for i = 1:size(LoadNames_RSMs,1) %Identify old file description
    try
        [RSMsNameInfo] = fMRIPrep_Name2Info(LoadNames_RSMs{i,1}{1,1});
        oldDesc=[oldDesc,RSMsNameInfo.desc];
    catch
        continue
    end
    break
end
newDesc='desc-SplitVarRCA'; %set file description name
[SavePaths_RCA] = strrepCell(filePaths_RSMs,'/SplitVarRSMs/',['/',AnalysisType,'/']);
[SavePaths_RCA] = strrepCell(SavePaths_RCA,['/',RSMsName,'/'],['/',AnalysisName,'/']);

[SaveDirs_RCA,SaveFileNames_RCA] = filepath2DirAndFileName(SavePaths_RCA);
[SaveFileNames_RCA] = strrepCell(SaveFileNames_RCA,oldDesc,newDesc);

for parcelNum=1:numParcels
    tic
    parcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
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
    UseLabels=UseLabels(:);
    RSMs=[];
    UseInd=[];
    All_RF = cell(1,2);
    All_RC = cell(1,2);
    [RSMs1,~,~] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','Split1_RSMs',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','SplitVarRSMs',...
        'AnalysisName',RSMsName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',parcelName);
    if isempty(RSMs1)
        continue
    end
    RSMs1=squeeze(real(RSMs1));
    [RSMs2,~,UseInd] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','Split2_RSMs',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','SplitVarRSMs',...
        'AnalysisName',RSMsName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',parcelName);
    if isempty(RSMs2)
        continue
    end
    RSMs2=squeeze(real(RSMs2));
    try
        if ndims(RSMs2) == 3
        [ OutVars1,OutVars2,OutVarsDiff ] = FullRCA_SplitVar(RSMs1,RSMs2,'vertIn',1,...
            'DownSample',DownSample,...
            'DownSampleN',DownSampleN,...
            'DownSampleReps',DownSampleReps,...
            'RunPermute',RunPermute,...
            'PermuteReps',PermuteReps,...
            'OutputPermDists',OutputPermDists,...
            'zNorm',zNormRSM);
        else
        [ OutVars1,OutVars2,OutVarsDiff ] = FullRCA_SplitVar(RSMs1,RSMs2,'vertIn',0,...
            'DownSample',DownSample,...
            'DownSampleN',DownSampleN,...
            'DownSampleReps',DownSampleReps,...
            'RunPermute',RunPermute,...
            'PermuteReps',PermuteReps,...
            'OutputPermDists',OutputPermDists,...
            'zNorm',zNormRSM);

        end
    catch
        continue
    end
    All_RF{1,1}=single(OutVars1.RF);
    All_RC{1,1}=single(OutVars1.RCMat); 
    All_RF{1,2}=single(OutVars2.RF);
    All_RC{1,2}=single(OutVars2.RCMat);     
    
    [AllSplitData] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','SplitData',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','SplitVarRSMs',...
        'AnalysisName',RSMsName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',parcelName); 
    
    All_RFDiff=single(OutVarsDiff.RF);
    All_RCDiff=single(OutVarsDiff.RCMat);
    UseInd=find(UseInd);
    UseInd=UseInd(:);
    for loadInd = 1:length(UseInd)
        loadNum=UseInd(loadInd,1);        
        if zNormRCA == 0
            Split1_RF=All_RF{1,1}(:,loadInd);
            Split1_RF=array2table(Split1_RF','VariableNames',UseLabels);
            Split1_RC=All_RC{1,1}(:,:,loadInd);
            Split1_RC=array2table(Split1_RC,'VariableNames',UseLabels,'RowNames',UseLabels);

            Split2_RF=All_RF{1,2}(:,loadInd);
            Split2_RF=array2table(Split2_RF','VariableNames',UseLabels);
            Split2_RC=All_RC{1,2}(:,:,loadInd);
            Split2_RC=array2table(Split2_RC,'VariableNames',UseLabels,'RowNames',UseLabels);

            RFDiff=All_RFDiff(:,loadInd);
            RFDiff=array2table(RFDiff','VariableNames',UseLabels);
            RCDiff=All_RCDiff(:,:,loadInd);
            RCDiff=array2table(RCDiff,'VariableNames',UseLabels,'RowNames',UseLabels);
        else
            Split1_RF=All_RF{1,1}(:,loadInd);
            Split1_RF=nan_zscore(Split1_RF);
            All_RF{1,1}(:,loadInd)=Split1_RF;
            Split1_RF=array2table(Split1_RF','VariableNames',UseLabels);
            
            Split1_RC=mat2uppertriuvectormat(All_RC{1,1}(:,:,loadInd));
            Split1_RC=vertRSM2SymRSM(nan_zscore(Split1_RC),[],1);
            All_RC{1,1}(:,:,loadInd)=Split1_RC;
            Split1_RC=array2table(Split1_RC,'VariableNames',UseLabels,'RowNames',UseLabels);

            Split2_RF=All_RF{1,2}(:,loadInd);
            Split2_RF=nan_zscore(Split2_RF);  
            All_RF{1,2}(:,loadInd)=Split2_RF;
            Split2_RF=array2table(Split2_RF','VariableNames',UseLabels);
            
            Split2_RC=mat2uppertriuvectormat(All_RC{1,2}(:,:,loadInd));
            Split2_RC=vertRSM2SymRSM(nan_zscore(Split2_RC),[],1); 
            All_RC{1,2}(:,:,loadInd)=Split2_RC;
            Split2_RC=array2table(Split2_RC,'VariableNames',UseLabels,'RowNames',UseLabels);
            
            RFDiff=table2array(Split2_RF)-table2array(Split1_RF);
            All_RFDiff(:,loadInd)=RFDiff;
            RFDiff=array2table(RFDiff,'VariableNames',UseLabels);
            
            RCDiff=table2array(Split2_RC)-table2array(Split1_RC);
            All_RCDiff(:,:,loadInd)=RCDiff;
            RCDiff=array2table(RCDiff,'VariableNames',UseLabels,'RowNames',UseLabels);
        end
        if ~exist(SaveDirs_RCA.(parcelName){loadNum,1},'file')
            mkdir(SaveDirs_RCA.(parcelName){loadNum,1});
        end
        SaveName=[SaveDirs_RCA.(parcelName){loadNum,1},SaveFileNames_RCA.(parcelName){loadNum,1}];
    
        % Split1 degree vars
        [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(Split1_RC),'ThresholdType','rank','Threshold',0.5);
        rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
        for nameNum=1:length(rowNames50)
            rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
        end
        tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);       
        [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(Split1_RC),'ThresholdType','rank','Threshold',0.25);
        rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
        for nameNum=1:length(rowNames25)
            rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
        end        
        tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);        
        [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(Split1_RC),'ThresholdType','rank','Threshold',0.10);
        rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
        for nameNum=1:length(rowNames10)
            rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
        end        
        tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);       
        [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(Split1_RC),'ThresholdType','rank','Threshold',0.05);
        rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
        for nameNum=1:length(rowNames05)
            rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
        end        
        tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);        
        Split1_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
        Split1_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
        Split1_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
        
        % Split2 degree vars
        [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(Split2_RC),'ThresholdType','rank','Threshold',0.5);
        rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
        for nameNum=1:length(rowNames50)
            rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
        end
        tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);       
        [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(Split2_RC),'ThresholdType','rank','Threshold',0.25);
        rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
        for nameNum=1:length(rowNames25)
            rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
        end        
        tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);        
        [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(Split2_RC),'ThresholdType','rank','Threshold',0.10);
        rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
        for nameNum=1:length(rowNames10)
            rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
        end        
        tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);       
        [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(Split2_RC),'ThresholdType','rank','Threshold',0.05);
        rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
        for nameNum=1:length(rowNames05)
            rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
        end        
        tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);        
        Split2_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
        Split2_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
        Split2_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
        
         % SplitDiff degree vars
        [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(RCDiff),'ThresholdType','rank','Threshold',0.5);
        rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
        for nameNum=1:length(rowNames50)
            rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
        end
        tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);       
        [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(RCDiff),'ThresholdType','rank','Threshold',0.25);
        rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
        for nameNum=1:length(rowNames25)
            rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
        end        
        tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);        
        [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(RCDiff),'ThresholdType','rank','Threshold',0.10);
        rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
        for nameNum=1:length(rowNames10)
            rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
        end        
        tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);       
        [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(RCDiff),'ThresholdType','rank','Threshold',0.05);
        rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
        for nameNum=1:length(rowNames05)
            rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
        end        
        tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);        
        SplitDiff_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
        SplitDiff_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
        SplitDiff_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
        
        SplitData=AllSplitData(loadInd,:);
        save(SaveName,'Split1_RF','Split2_RF','Split1_RC','Split2_RC','RFDiff','RCDiff','Split1_RC_Degree','Split2_RC_Degree','SplitDiff_RC_Degree','Split1_RC_WholeCM','Split2_RC_WholeCM','SplitDiff_RC_WholeCM','SplitData','AnalysisParameters');
    end

    [Split1_Group_RF,Split1_Group_RC,Split1_Group_RCMat_Ts,Split1_Group_RCMat_Means,Split1_RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(All_RF{1,1},All_RC{1,1},GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
        'vecName','Split1_RF','matName','Split1_RC','SavePrefix','Split1_','permZMat',OutVars1.RCMat_permZ,'permPMat',OutVars1.RCMat_permP,'permZVec',OutVars1.RF_permZ,'permPVec',OutVars1.RF_permP,'permZMat_Corrected',OutVars1.RCMat_permZ_Corrected,'permPMat_Corrected',OutVars1.RCMat_permP_Corrected,'permZVec_Corrected',OutVars1.RF_permZ_Corrected,'permPVec_Corrected',OutVars1.RF_permP_Corrected); 
    if ComputeReliability==1
        Split1_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
        Split1_RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
        Split1_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
        Split1_RC_losoReliability=ReliabilityResults.Mat_losoReliability;
    end
    [Split2_Group_RF,Split2_Group_RC,Split2_Group_RCMat_Ts,Split2_Group_RCMat_Means,Split2_RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(All_RF{1,2},All_RC{1,2},GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
       'vecName','Split2_RF','matName','Split2_RC','SavePrefix','Split2_','permZMat',OutVars2.RCMat_permZ,'permPMat',OutVars2.RCMat_permP,'permZVec',OutVars2.RF_permZ,'permPVec',OutVars2.RF_permP,'permZMat_Corrected',OutVars2.RCMat_permZ_Corrected,'permPMat_Corrected',OutVars2.RCMat_permP_Corrected,'permZVec_Corrected',OutVars2.RF_permZ_Corrected,'permPVec_Corrected',OutVars2.RF_permP_Corrected);
    if ComputeReliability==1
        Split2_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
        Split2_RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
        Split2_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
        Split2_RC_losoReliability=ReliabilityResults.Mat_losoReliability;    
    end
    [SplitDiff_Group_RF,SplitDiff_Group_RC,SplitDiff_Group_RCMat_Ts,SplitDiff_Group_RCMat_Means,SplitDiff_RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(All_RFDiff,All_RCDiff,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
       'vecName','SplitDiff_RF','matName','SplitDiff_RC','SavePrefix','SplitDiff_','permZMat',OutVarsDiff.RCMat_permZ,'permPMat',OutVarsDiff.RCMat_permP,'permZVec',OutVarsDiff.RF_permZ,'permPVec',OutVarsDiff.RF_permP,'permZMat_Corrected',OutVarsDiff.RCMat_permZ_Corrected,'permPMat_Corrected',OutVarsDiff.RCMat_permP_Corrected,'permZVec_Corrected',OutVarsDiff.RF_permZ_Corrected,'permPVec_Corrected',OutVarsDiff.RF_permP_Corrected);
    if ComputeReliability==1
        SplitDiff_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
        SplitDiff_RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
        SplitDiff_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
        SplitDiff_RC_losoReliability=ReliabilityResults.Mat_losoReliability;  
    end
    SplitDiffs=mean(AllSplitData,1);
    if length(UseLabels)<350
        %% compute parcellation and SS-based RSMs for 3rd order RSA.
        RFMats=All_RF{1,1};
        RCMats=All_RC{1,1};
        RCMatsVert=mat2uppertriuvectormat(RCMats);
        ThirdOrderRSMs.Split1.ssRSM_RF=corrcoef(RFMats);
        ThirdOrderRSMs.Split1.parcelRSM_RF=corrcoef(RFMats');
        ThirdOrderRSMs.Split1.ssRSM_RC=corrcoef(RCMatsVert);
        ThirdOrderRSMs.Split1.parcelRSM_RC=corrcoef(RCMatsVert');
        ThirdOrderRSMs.Split1.ssRSM_RF_scaled=corrcoef(scaleVals(RFMats,2));
        ThirdOrderRSMs.Split1.parcelRSM_RF_scaled=corrcoef(scaleVals(RFMats',2));  
        ThirdOrderRSMs.Split1.ssRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert,2));
        ThirdOrderRSMs.Split1.parcelRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert',2));  

        RFMats=All_RF{1,2};
        RCMats=All_RC{1,2};
        RCMatsVert=mat2uppertriuvectormat(RCMats);
        ThirdOrderRSMs.Split2.ssRSM_RF=corrcoef(RFMats);
        ThirdOrderRSMs.Split2.parcelRSM_RF=corrcoef(RFMats');
        ThirdOrderRSMs.Split2.ssRSM_RC=corrcoef(RCMatsVert);
        ThirdOrderRSMs.Split2.parcelRSM_RC=corrcoef(RCMatsVert');
        ThirdOrderRSMs.Split2.ssRSM_RF_scaled=corrcoef(scaleVals(RFMats,2));
        ThirdOrderRSMs.Split2.parcelRSM_RF_scaled=corrcoef(scaleVals(RFMats',2));  
        ThirdOrderRSMs.Split2.ssRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert,2));
        ThirdOrderRSMs.Split2.parcelRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert',2)); 

        RFMats=All_RFDiff;
        RCMats=All_RCDiff;
        RCMatsVert=mat2uppertriuvectormat(RCMats);
        ThirdOrderRSMs.SplitDiff.ssRSM_RF=corrcoef(RFMats);
        ThirdOrderRSMs.SplitDiff.parcelRSM_RF=corrcoef(RFMats');
        ThirdOrderRSMs.SplitDiff.ssRSM_RC=corrcoef(RCMatsVert);
        ThirdOrderRSMs.SplitDiff.parcelRSM_RC=corrcoef(RCMatsVert');
        ThirdOrderRSMs.SplitDiff.ssRSM_RF_scaled=corrcoef(scaleVals(RFMats,2));
        ThirdOrderRSMs.SplitDiff.parcelRSM_RF_scaled=corrcoef(scaleVals(RFMats',2));  
        ThirdOrderRSMs.SplitDiff.ssRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert,2));
        ThirdOrderRSMs.SplitDiff.parcelRSM_RC_scaled=corrcoef(scaleVals(RCMatsVert',2));         
        UseTable=fmriprep_table(UseInd,:);
        save([GroupAnalysisDir,'ThirdOrderRSMs'],'ThirdOrderRSMs','UseInd','UseTable');
    end
    save([GroupAnalysisDir,'AnalysisParams_RCA'],'AnalysisParameters');
    save([GroupAnalysisDir,'Group_RF'],'Split1_Group_RF','Split2_Group_RF','SplitDiff_Group_RF','SplitDiffs');
    save([GroupAnalysisDir,'Group_RC'],'Split1_Group_RC','Split2_Group_RC','SplitDiff_Group_RC','SplitDiffs');
    save([GroupAnalysisDir,'Group_RCMat'],'Split1_Group_RCMat_Ts','Split1_Group_RCMat_Means','Split1_RC_Degree_Table','Split2_Group_RCMat_Ts','Split2_Group_RCMat_Means','Split2_RC_Degree_Table','SplitDiff_Group_RCMat_Ts','SplitDiff_Group_RCMat_Means','SplitDiff_RC_Degree_Table');  
    if ComputeReliability==1
        save([GroupAnalysisDir,'Group_ParcelReliability'],'Split1_RF_SplitHalfReliability','Split1_RC_SplitHalfReliability','Split1_RF_losoReliability','Split1_RC_losoReliability','Split2_RF_SplitHalfReliability','Split2_RC_SplitHalfReliability','Split2_RF_losoReliability','Split2_RC_losoReliability','SplitDiff_RF_SplitHalfReliability','SplitDiff_RC_SplitHalfReliability','SplitDiff_RF_losoReliability','SplitDiff_RC_losoReliability');  
    end
    toc
end
end

end

