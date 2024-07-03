function [AnalysisParameters] = ComputeSearchlight_SeedRC(fmriprep_table,ExperimentsDir,varargin)
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
[SeedRSMType] = VariableSetter('SeedRSMType',[],varargin);
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
CleanRSMs = VariableSetter( 'CleanRSMs',[],varargin);
RefInfo = VariableSetter( 'RefInfo',[],varargin);
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
    SeedRSMType=AnalysisParameters.SeedRSMType;
    zNormRSM=AnalysisParameters.zNormRSM;
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
AnalysisType='SeedRC';
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
if isempty(SeedRSMType)
    SingleSelect=1;
    [SeedRSMType] = uiNameSelect({'Parcel','Reference'},'Select type of seed RSM:',SingleSelect);
end
if strcmpi(SeedRSMType,'Reference')
    if isempty(CleanRSMs)
        CleanRSMs=uiNameSelect({'Yes','No'},'Use cleaned RSMs?');
    end
    if strcmpi(CleanRSMs,'Yes')
        clAppend='_cl';
    else
        clAppend='';
    end   
    if isempty(RefInfo)
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','GroupType','GroupAnalysis','FileNames','ParcelRSMs.mat');
    else
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','FileNames','ParcelRSMs.mat','GroupType','GroupAnalysis',...
            'AnalysisNames',RefInfo.AnalysisNames,'TableNames',RefInfo.TableNames,'ParcellationNames',RefInfo.ParcellationNames,'LoadVarNames_Parcel',RefInfo.LoadVarNames_Parcel);
    end
    RefParcelNames=RefInfo.ParcellationNames;
    if strcmpi(CleanRSMs,'Yes')
        AllRSMs=load(AllPaths_RefRSMs{1,1},'ParcelRSMs_p50','ParcelZRSMs_p50');
        AllRSMs.ParcelRSMs=AllRSMs.ParcelRSMs_p50;
        AllRSMs.ParcelZRSMs=AllRSMs.ParcelZRSMs_p50;
    else
        AllRSMs=load(AllPaths_RefRSMs{1,1},'ParcelRSMs','ParcelZRSMs');
    end  
    SingleSelect=1;
    SeedName=uiNameSelect(AllRSMs.ParcelRSMs.Properties.VariableNames,'Select Parcel ROI for Seed Analysis',SingleSelect);
    SeedRSMs=AllRSMs.ParcelRSMs.(SeedName);
    ParcelName=RefInfo.ParcellationNames{1,1};
    parcelName=ParcelName;    
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
    SeedMask=single(UseMask==find(ismember(UseLabels,SeedName)));    
    RefAppend=[RefInfo.TableNames{1,1},'_',RefInfo.AnalysisNames{1,1}];
    RefAppend=strrep(RefAppend,'_RSMs','');
    SeedAppend=[RefAppend,clAppend];
    SeedName=[SearchlightName,'_',ParcelName,'_',SeedName];
else
    [SeedRSMs,fmriprep_table,selectind,DataLabels,~,~,~,LoadParams] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','RSMs',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','RSMs',...
        'AnalysisName',[],...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',[]);
    SeedName=uiNameSelect(DataLabels,'Select Parcel ROI for Seed Analysis',1);
    SeedRSMs=squeeze(SeedRSMs(:,find(ismember(DataLabels,SeedName)),:));
    SeedRSMs=repmat(nansum(SeedRSMs,2)/(size(SeedRSMs,2)-1),[1,size(SeedRSMs,2)])-(SeedRSMs/(size(SeedRSMs,2)-1));
    SeedAppend=LoadParams.AnalysisName;
    ParcelName=LoadParams.ParcelName;

    load(['Parcellations/',ParcelName],'UseLabels','UseMask');
    UseLabels=UseLabels(:);
    SeedMask=single(UseMask==find(ismember(UseLabels,SeedName)));
    SeedName=[SearchlightName,'_',ParcelName,'_',SeedName];  
end
SeedName=strrep(SeedName,'Searchlight_Shape-','');
AnalysisParameters.SeedRSMType=SeedRSMType;
if isempty(AnalysisName)
    RCAName=strrep(RSMsName,'RSMs_',['SeedRC_',fmriprep_table_name,'_']);
    if DefaultName~=1
        AnalysisName=uiEnterName([RCAName,'_',SeedAppend],['Enter name for ',AnalysisType,newline,'analysis below:']);
    else
       AnalysisName=[RCAName,'_',SeedAppend];
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
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SeedName,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SeedName,'/'];
else
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SeedName,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SeedName,'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);     
end
if ~exist(GroupBrainMapsDir,'file')
    mkdir(GroupBrainMapsDir);     
end   

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
    SaveDir=[fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/',SeedName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end
    descript1='desc-SeedRC'; %set file description name
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
    if strcmpi(SeedRSMType,'Parcel')
        sRSM=SeedRSMs(:,dataInd);
    else
        sRSM=SeedRSMs;
    end
    slRSMs=load(filePaths_Searchlight{dataInd,1});
    if ~any(isnan(slRSMs.RSMs(:)))
        seedData=corr(slRSMs.RSMs,sRSM);
    else
        seedData=corr(slRSMs.RSMs,sRSM,'rows','pairwise');
    end
    seedMap=single(slRSMs.rsm_mask);
    useInd=find(seedMap);
    seedMap(seedMap==0)=nan;
    seedMap(useInd)=seedData;
    save(SaveNames{1,1},'seedMap','AnalysisParameters');    
end
[seedMaps] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','seedMap',...
    'LoadVarFormat','Array',...
    'TableOrArray','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','SeedRC',...
    'AnalysisName',AnalysisName,...
    'SubjectOrRun',SubjectOrRun,...
    'ParcelName',SeedName);
ThreshMask=single(sum(single(~isnan(seedMaps)),4)>round(0.75*size(seedMaps,4)));
ThreshMask(SeedMask==1)=0;
ThreshMask=repmat(ThreshMask,[1,1,1,size(seedMaps,4)]);
seedMaps(ThreshMask==0)=nan;
[groupSeedRC_t,groupSeedRC_z,groupSeedRC_mean]=getTval(atanh(seedMaps),4);
groupSeedRC_t(isinf(groupSeedRC_t))=nan;
groupSeedRC_z(isinf(groupSeedRC_z))=nan;
groupSeedRC_mean(isinf(groupSeedRC_mean))=nan;
BrainMap=cat(4,groupSeedRC_mean,groupSeedRC_t,groupSeedRC_z);
SaveBrik_3mmMNI(BrainMap,{'group_mean','group_t','group_z'},[GroupBrainMapsDir,SeedName,'_SearchlightMaps']);
save([GroupAnalysisDir,SeedName,'_SearchlightMaps'],'groupSeedRC_mean','groupSeedRC_t','groupSeedRC_z','AnalysisParameters');      
%%Run ClusterThresholds
% if strcmpi(ComputeClusterThresh,'Indv') || strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
%     if strcmpi(ComputeClusterThresh,'Indv') || strcmpi(ComputeClusterThresh,'Both')
%         IndvClusterDists=cell(size(AllRSMs,3),1);
%     end
%     if strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
%         GroupClusterDists=nan(ClusterThreshReps,3);
%     end    
%     for rep = 1:ClusterThreshReps
%         disp(['Running cluster correction: ',num2str((rep/ClusterThreshReps)*100),'% Complete.'])
%         tempAllRF=cell(size(totalMask,1),1);
%         for i = 1:size(AllRSMs,3)
%             AllRSMs(:,:,i)=modelscramble(AllRSMs(:,:,i),'VH_once',[],'vertIn',1);
%         end
%         parfor i = 1:size(AllRSMs,2)    
%             tempRSMs=AllRSMs(:,i,:);
%             nanInd=isnan(tempRSMs(1,1,:));
%             tempRF=nan(1,length(nanInd));
%             [ OutVars ] = FullRCA( tempRSMs(:,1,nanInd==0),'vertIn',1,...
%                 'DownSample',DownSample,...
%                 'DownSampleN',DownSampleN,...
%                 'DownSampleReps',DownSampleReps,...
%                 'RunPermute',RunPermute,...
%                 'PermuteReps',PermuteReps,...
%                 'OutputPermDists',OutputPermDists,...
%                 'zNorm',zNormRSM);    
%             tempRF(1,nanInd==0)=single(OutVars.RF);
%             tempAllRF{i,1}=tempRF;
%         end  
%         tempAllRF=single(squeeze(cell2nDMAT(tempAllRF))); 
%         if strcmpi(ComputeClusterThresh,'Group') || strcmpi(ComputeClusterThresh,'Both')
%             [tempgroupRF_t,tempgroupRF_z]=getTval(atanh(tempAllRF),1);        
%             tempMask=single(nan(length(groupMaskInd),1));
%             tempMask(groupMaskInd)=tempgroupRF_t;       
%             tempgroupRF_t=tempMask;        
%             tempMask=single(nan(length(groupMaskInd),1));
%             tempMask(groupMaskInd)=tempgroupRF_z;
%             tempgroupRF_z=tempMask;
%             tempgroupRF_z(isinf(tempgroupRF_z))=tempgroupRF_t(isinf(tempgroupRF_z));
%             tempgroupRF_z=reshape(tempgroupRF_z,mapSize);
%             [maskcoords,vals]=mat2coords(tempgroupRF_z);        
%             [ Clusters ] = ClusterFinder(vals,maskcoords,'pos','stat',UncorrectedThresholds,1,'side');
%             GroupClusterDists(rep,:)=Clusters.MaxClusterSizes;
%         end           
%     end
%     GroupClusterThreshs=prctile(GroupClusterDists,95);    
%     [maskcoords,vals]=mat2coords(groupRF_z);  
%     ThreshGroupRF_z=[];
%     mapLabels=cell(1,length(UncorrectedThresholds));
%     Clusters=cell(length(UncorrectedThresholds),1);
%     for j = 1:length(UncorrectedThresholds)
%         [ Clusters{j,1} ] = ClusterFinder(vals,maskcoords,'pos','stat',UncorrectedThresholds(1,j),GroupClusterThreshs(1,j),'side');
%         threshmat=coords2mat(Clusters{j,1}.AllThresholdedCoords{1,1},groupRF_z,ones(length(Clusters{j,1}.AllThresholdedCoords{1,1}),1));
%         mapLabels{1,j}=['VoxThresh_Zgt',num2str4filename(UncorrectedThresholds(1,j),2)];
%         ThreshGroupRF_z=cat(4,ThreshGroupRF_z,groupRF_z.*threshmat);
%     end
%     SaveBrik_3mmMNI(ThreshGroupRF_z,mapLabels,[GroupBrainMapsDir,'groupRF_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)]);
%     save([GroupAnalysisDir,'groupRF_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)],'ThreshGroupRF_z','GroupClusterDists','GroupClusterThreshs','UncorrectedThresholds','Clusters','AnalysisParameters');      
% end


