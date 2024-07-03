function [AnalysisParameters] = ComputeSearchlightRCA(fmriprep_table,ExperimentsDir,varargin)
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
[AnalysisName] = VariableSetter('AnalysisName',[''],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[RSMsName] = VariableSetter('RSMsName',[],varargin);
[SearchlightName] = VariableSetter('SearchlightName',[],varargin);
[nThreshold] = VariableSetter('nThreshold',[],varargin);
[ComputeClusterThresh] = VariableSetter('ComputeClusterThresh',[],varargin);
[UncorrectedThresholds] = VariableSetter('UncorrectedThresholds',[2.576,2.807,3.2905],varargin);
[ClusterThreshReps] = VariableSetter('ClusterThreshReps',1000,varargin);

%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
[DownSample] = VariableSetter('DownSample',[],varargin);
%Subject or run level analysis. Will prompt request.
[DownSampleN] = VariableSetter('DownSampleN',[],varargin);
%Subject or run level analysis. Will prompt request.
[DownSampleReps] = VariableSetter('DownSampleReps',100,varargin);
%Subject or run level analysis. Will prompt request.
[RunPermute] = VariableSetter('RunPermute',[],varargin);
%Subject or run level analysis. Will prompt request.
[PermuteReps] = VariableSetter('PermuteReps',[],varargin);
%Subject or run level analysis. Will prompt request.
[OutputPermDists] = VariableSetter('OutputPermDists',0,varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
zNormRSM = VariableSetter( 'zNormRSM',[],varargin);
RSMGroup = VariableSetter( 'RSMGroup',[],varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_RCA,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','RCA','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_RCA.Properties.VariableNames;
        catch
            filePaths_RCA=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_RCA)
            for j = 1:height(filePaths_RCA)
                try
                    load(filePaths_RCA.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
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
    zNormRSM=AnalysisParameters.zNormRSM;
    RSMGroup=AnalysisParameters.RSMGroup;
    DownSample=AnalysisParameters.DownSample;
    DownSampleN=AnalysisParameters.DownSampleN;
    DownSampleReps=AnalysisParameters.DownSampleReps;
    RunPermute=AnalysisParameters.RunPermute;
    RSMsName=AnalysisParameters.RSMsName;
    PermuteReps=AnalysisParameters.PermuteReps;
    OutputPermDists=AnalysisParameters.OutputPermDists;
    SearchlightName=AnalysisParameters.SearchlightName;
    ComputeClusterThresh=AnalysisParameters.ComputeClusterThresh;
    fmriprep_table_name=AnalysisParameters.fmriprep_table_name;  
else
    AnalysisParameters=struct;
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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='RCA';
%Allows you to set name for this particular analysis
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
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
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;

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
    zNormSuffix='_zNorm';
else
    zNormSuffix='';
end
AnalysisParameters.zNormRSM=zNormRSM;
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
RSMGroupSuffix='';
if isempty(RSMGroup)
    SingleSelect=1; %Allows only a single value to be selected.
    [RSMGroup] = uiNameSelect({'All','W','X'},'Select section of RSM to use:',SingleSelect);
    if ~strcmpi(RSMGroup,'All')
        RSMGroupSuffix=['_',RSMGroup];
    end
end
%% Set searchligh setting
if isempty(SearchlightName)
    [filePaths_Searchlight,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',RSMsName,'TitleTextName','Select fMRI input for GLM:');
    SearchlightName=uiNameSelect(filePaths_Searchlight.Properties.VariableNames,'Select searchlight data:');
    if iscell(SearchlightName)
        SearchlightName=SearchlightName{1,1};
    end
end
AnalysisParameters.SearchlightName=SearchlightName;
if isempty(ComputeClusterThresh)
    SingleSelect=1; %Allows only a single value to be selected.
    [ComputeClusterThresh] = uiNameSelect({'None','Indv','Group','Both'},'Compute label scramble cluster size correction:',SingleSelect);
end
AnalysisParameters.ComputeClusterThresh=ComputeClusterThresh;
[totalMask,fmriprep_table,selectind] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','rsm_mask',...
    'LoadVarFormat','Array',...
    'TableOrArray','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','RSMs',...
    'AnalysisName',RSMsName,...
    'SubjectOrRun',SubjectOrRun,...
    'ParcelName',SearchlightName);

if isempty(AnalysisName)
    RCAName=strrep(RSMsName,'RSMs_',['RCA_',fmriprep_table_name,'_']);
    if DefaultName~=1
        AnalysisName=uiEnterName([RCAName,zNormSuffix,DownSampleSuffix,PermuteSuffix,RSMGroupSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    else
       AnalysisName=[RCAName,DownSampleSuffix,PermuteSuffix];
    end
end
[filePaths_Searchlight] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',RSMsName,'TitleTextName','Select fMRI input for GLM:');
filePaths_Searchlight=filePaths_Searchlight.(SearchlightName);
AnalysisParameters.SubjectOrRun=SubjectOrRun;
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
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
if bySS == 1
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightName,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightName,'/'];
else
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SearchlightName,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SearchlightName,'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);     
end
if ~exist(GroupBrainMapsDir,'file')
    mkdir(GroupBrainMapsDir);     
end   
groupMask=sum(totalMask,4);
mapSize=size(groupMask);
groupMask=groupMask(:);
for i = 1:size(totalMask,4)
    tempMask=totalMask(:,:,:,i);
    tempMask(tempMask==1)=1:sum(tempMask(:));
    totalMask(:,:,:,i)=tempMask;
end
totalMask=reshape(totalMask,[size(totalMask,1)*size(totalMask,2)*size(totalMask,3),size(totalMask,4)]);
if isempty(nThreshold)
    nThreshold=round(sum(single(selectind))*0.5);
end
groupMaskInd=groupMask>nThreshold;
totalMask=totalMask(groupMaskInd,:);
groupMask=groupMask(groupMaskInd,:);
AllRSMs=[];
for i =1:size(filePaths_Searchlight,1)
    loadData=load(filePaths_Searchlight{i,1},'RSMs');
    tempRSMs=loadData.RSMs;
    if isempty(AllRSMs)
        AllRSMs=nan(size(tempRSMs,1),size(totalMask,2),size(totalMask,1),'single');
    end
    tempRSMs=[tempRSMs,nan(size(tempRSMs,1),1)];
    selectind=totalMask(:,i);
    selectind(selectind==0,:)=size(tempRSMs,2);
    AllRSMs(:,i,:)=tempRSMs(:,selectind);
end
AllRSMs=permute(AllRSMs,[1,3,2]);
%AllRF=nan(size(totalMask,1),size(totalMask,2),'single');
AllRF=cell(size(totalMask,1),1);
if strcmpi(RSMGroup,'W')
    [useRSMs,~] = RSM2RSMGroups(AllRSMs,1);
elseif strcmpi(RSMGroup,'X')
    [~,useRSMs] = RSM2RSMGroups(AllRSMs,1);
else
    useRSMs=AllRSMs;
end
tic
parfor i = 1:size(useRSMs,2)    
    tempRSMs=useRSMs(:,i,:);
    nanInd=isnan(tempRSMs(1,1,:));
    tempRF=nan(1,length(nanInd));
    [ OutVars ] = FullRCA( tempRSMs(:,1,nanInd==0),'vertIn',1,...
        'DownSample',DownSample,...
        'DownSampleN',DownSampleN,...
        'DownSampleReps',DownSampleReps,...
        'RunPermute',RunPermute,...
        'PermuteReps',PermuteReps,...
        'OutputPermDists',OutputPermDists,...
        'zNorm',zNormRSM,...
        'RSMGroup','All');    
    tempRF(1,nanInd==0)=single(OutVars.RF);
    AllRF{i,1}=tempRF;
end    
toc
AllRF=single(squeeze(cell2nDMAT(AllRF)));
[groupRF_t,groupRF_z,groupRF_mean]=getTval(atanh(AllRF),1);
tempMask=single(nan(length(groupMaskInd),1));
tempMask(groupMaskInd)=groupRF_t;
groupRF_t=tempMask;
tempMask=single(nan(length(groupMaskInd),1));
tempMask(groupMaskInd)=groupRF_z;
groupRF_z=tempMask;
groupRF_z(isinf(groupRF_z))=groupRF_t(isinf(groupRF_z));
tempMask=single(nan(length(groupMaskInd),1));
tempMask(groupMaskInd)=groupRF_mean;
groupRF_mean=tempMask;
groupRF_t=reshape(groupRF_t,mapSize);
groupRF_z=reshape(groupRF_z,mapSize);
groupRF_mean=reshape(groupRF_mean,mapSize);

BrainMap=cat(4,groupRF_mean,groupRF_t,groupRF_z);
SaveBrik_3mmMNI(BrainMap,{'groupRF_mean','groupRF_t','groupRF_z'},[GroupBrainMapsDir,'groupRF_SearchlightMaps']);
save([GroupAnalysisDir,'groupRF_SearchlightMaps'],'groupRF_mean','groupRF_t','groupRF_z','AnalysisParameters');      

%%Save SSMaps
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
    SaveDir=[fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/',SearchlightName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end
    descript1='desc-RCA'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(1);
    SavePrefix=[ExperimentsDir,SaveDir,'/'];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    else  
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
            disp(['Skipping-- all files exist: ',SaveName]);
            continue
        end
    end  
    RFMap=single(nan(length(groupMaskInd),1));
    RFMap(groupMaskInd)=AllRF(dataInd,:);
    RFMap=reshape(RFMap,mapSize);    
 
    save(SaveNames{1,1},'RFMap','AnalysisParameters');    
end

%%Run ClusterThresholds
if strcmpi(ComputeClusterThresh,'Indv') || strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
    if strcmpi(ComputeClusterThresh,'Indv') || strcmpi(ComputeClusterThresh,'Both')
        IndvClusterDists=cell(size(AllRSMs,3),1);
    end
    if strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
        GroupClusterDists=nan(ClusterThreshReps,3);
    end    
    for rep = 1:ClusterThreshReps
        disp(['Running cluster correction: ',num2str((rep/ClusterThreshReps)*100),'% Complete.'])
        tempAllRF=cell(size(totalMask,1),1);
        for i = 1:size(AllRSMs,3)
            AllRSMs(:,:,i)=modelscramble(AllRSMs(:,:,i),'VH_once',[],'vertIn',1);
        end
        if strcmpi(RSMGroup,'W')
            [useRSMs,~] = RSM2RSMGroups(AllRSMs,1);
        elseif strcmpi(RSMGroup,'X')
            [~,useRSMs] = RSM2RSMGroups(AllRSMs,1);
        end        
        parfor i = 1:size(useRSMs,2)    
            tempRSMs=useRSMs(:,i,:);
            nanInd=isnan(tempRSMs(1,1,:));
            tempRF=nan(1,length(nanInd));
            [ OutVars ] = FullRCA( tempRSMs(:,1,nanInd==0),'vertIn',1,...
                'DownSample',DownSample,...
                'DownSampleN',DownSampleN,...
                'DownSampleReps',DownSampleReps,...
                'RunPermute',RunPermute,...
                'PermuteReps',PermuteReps,...
                'OutputPermDists',OutputPermDists,...
                'zNorm',zNormRSM,...
                'RSMGroup','All');    
            tempRF(1,nanInd==0)=single(OutVars.RF);
            tempAllRF{i,1}=tempRF;
        end  
        tempAllRF=single(squeeze(cell2nDMAT(tempAllRF))); 
        if strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
            [tempgroupRF_t,tempgroupRF_z]=getTval(atanh(tempAllRF),1);        
            tempMask=single(nan(length(groupMaskInd),1));
            tempMask(groupMaskInd)=tempgroupRF_t;       
            tempgroupRF_t=tempMask;        
            tempMask=single(nan(length(groupMaskInd),1));
            tempMask(groupMaskInd)=tempgroupRF_z;
            tempgroupRF_z=tempMask;
            tempgroupRF_z(isinf(tempgroupRF_z))=tempgroupRF_t(isinf(tempgroupRF_z));
            tempgroupRF_z=reshape(tempgroupRF_z,mapSize);
            [maskcoords,vals]=mat2coords(tempgroupRF_z);        
            [ Clusters ] = ClusterFinder(vals,maskcoords,'pos','stat',UncorrectedThresholds,1,'side');
            GroupClusterDists(rep,:)=Clusters.MaxClusterSizes;
        end           
    end
    GroupClusterThreshs=prctile(GroupClusterDists,95);    
    [maskcoords,vals]=mat2coords(groupRF_z);  
    ThreshGroupRF_z=[];
    mapLabels=cell(1,length(UncorrectedThresholds));
    Clusters=cell(length(UncorrectedThresholds),1);
    for j = 1:length(UncorrectedThresholds)
        [ Clusters{j,1} ] = ClusterFinder(vals,maskcoords,'pos','stat',UncorrectedThresholds(1,j),GroupClusterThreshs(1,j),'side');
        threshmat=coords2mat(Clusters{j,1}.AllThresholdedCoords{1,1},groupRF_z,ones(length(Clusters{j,1}.AllThresholdedCoords{1,1}),1));
        mapLabels{1,j}=['VoxThresh_Zgt',num2str4filename(UncorrectedThresholds(1,j),2)];
        ThreshGroupRF_z=cat(4,ThreshGroupRF_z,groupRF_z.*threshmat);
    end
    SaveBrik_3mmMNI(ThreshGroupRF_z,mapLabels,[GroupBrainMapsDir,'groupRF_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)]);
    save([GroupAnalysisDir,'groupRF_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)],'ThreshGroupRF_z','GroupClusterDists','GroupClusterThreshs','UncorrectedThresholds','Clusters','AnalysisParameters');      
end


