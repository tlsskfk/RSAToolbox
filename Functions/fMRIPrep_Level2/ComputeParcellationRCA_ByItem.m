function [AnalysisParameters] = ComputeParcellationRCA_ByItem(fmriprep_table,ExperimentsDir,varargin)
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
%Subject or run level analysis. Will prompt request.
[RSMsName] = VariableSetter('RSMsName',[],varargin);
%Subject or run level analysis. Will prompt request.
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
%Subject or run level analysis. Will prompt request.
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
[OutputPermDists] = VariableSetter('OutputPermDists',[],varargin);
[ConditionNames] = VariableSetter('ConditionNames',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
zNormRSM = VariableSetter( 'zNormRSM',[],varargin);
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
    FigureInfo=AnalysisParameters.FigureInfo;
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    zNormRSM=AnalysisParameters.zNormRSM;
    DownSample=AnalysisParameters.DownSample;
    DownSampleN=AnalysisParameters.DownSampleN;
    DownSampleReps=AnalysisParameters.DownSampleReps;
    RunPermute=AnalysisParameters.RunPermute;
    RSMsName=AnalysisParameters.RSMsName;
    PermuteReps=AnalysisParameters.PermuteReps;
    OutputPermDists=AnalysisParameters.OutputPermDists;
    fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    ParcelNames=AnalysisParameters.ParcelNames;   
else
    AnalysisParameters=struct;
    ParcelNames=[];
end

if isempty(fmriprep_table_name)
[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
AnalysisParameters.fmriprep_table_name=fmriprep_table_name;
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='RCAbyItem';
%Allows you to set name for this particular analysis

%% Compile filepaths for input files for the analysis
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',RSMsName,'TitleTextName',['Select RSMs for',newline,'Representational Connectivity Analysis:']);
if isempty(AnalysisName)
    RCAName=strrep(RSMsName,'RSMs_',['RCAbyItem_',fmriprep_table_name,'_']);
    if DefaultName~=1
        AnalysisName=uiEnterName([RCAName,zNormSuffix,DownSampleSuffix,PermuteSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    else
       AnalysisName=[RCAName,DownSampleSuffix,PermuteSuffix];
    end
end

%% Identify parcellations availible and select ones to run analysis on.
if isempty(ParcelNames)
    ParcelNames=filePaths_RSMs.Properties.VariableNames;
    ParcelNames=ParcelNames(:);
    %ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
    ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
end
AnalysisParameters.ParcelNames=ParcelNames;
filePaths_RSMs=filePaths_RSMs(:,ParcelNames);
numParcels=length(ParcelNames);

for i = 1:height(filePaths_RSMs)
    tempLoad=filePaths_RSMs{i,1}{1,1};
    if ~isempty(tempLoad)
        break
    end
end
tempLoadVars=load(tempLoad,'AnalysisParameters');
if isempty(tempLoadVars.AnalysisParameters.ConditionNames)
    [filePaths_ActPat,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ActivationPatterns','AnalysisName',tempLoadVars.AnalysisParameters.ActivationPatternName,'TitleTextName',['Select RSMs for',newline,'Representational Connectivity Analysis:']);
    for i = 1:height(filePaths_ActPat)
        tempLoad=filePaths_ActPat{i,1}{1,1};
        if ~isempty(tempLoad)
            break
        end
    end
    tempLoadVars=load(tempLoad,'AnalysisParameters');
    ConditionNames=tempLoadVars.AnalysisParameters.ConditionNames;    
else
    ConditionNames=tempLoadVars.AnalysisParameters.ConditionNames;
end

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
newDesc='desc-RCA_byItem'; %set file description name
[SavePaths_RCA] = strrepCell(filePaths_RSMs,'/RSMs/',['/',AnalysisType,'/']);
[SavePaths_RCA] = strrepCell(SavePaths_RCA,['/',RSMsName,'/'],['/',AnalysisName,'/']);
[SavePaths_RCA] = strrepCell(SavePaths_RCA,['/',RSMsName,'/'],['/',AnalysisName,'/']);
[SaveDirs_RCA,SaveFileNames_RCA] = filepath2DirAndFileName(SavePaths_RCA);
[SaveFileNames_RCA] = strrepCell(SaveFileNames_RCA,oldDesc,newDesc);

if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end
AnalysisParameters.AnalysisType=AnalysisType;
AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.RSMsName=RSMsName;

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress


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
    AllRSMs=[];
    UseInd=[];
    for loadNum = 1:size(LoadNames_RSMs,1) 
        loadPath=filePaths_RSMs.(parcelName){loadNum,1};
        if isempty(loadPath)
            continue
        else
            try
                load(loadPath,'RSMs');
            catch
                continue
            end
        end
        try          
            AllRSMs=cat(4,AllRSMs,RSMs);
        catch
            continue
        end
        UseInd=[UseInd;loadNum];
    end
    if isempty(AllRSMs)
        continue
    end
    UseTable=fmriprep_table(UseInd,:);
    AllRSMs=squeeze(real(AllRSMs)); 
    if ndims(AllRSMs) == 2
        AllRSMs=reshape(AllRSMs,[size(AllRSMs,1),1,size(AllRSMs,2)]);
    end    
    AllRSMs=single(AllRSMs);
    RSMSize = vertRSM2SymRSM( ones(size(AllRSMs,1),1));
    RSMSize=size(RSMSize,1);
    SubSample=ones(RSMSize,1);
    SubRSMVertSize=((RSMSize-1)^2-(RSMSize-1))/2;
    RFMats=[];
    RCMats=[];
    ALLRF=zeros(size(AllRSMs,2),size(AllRSMs,3),RSMSize,'single');
    ALLRC=zeros(size(AllRSMs,2),size(AllRSMs,2),size(AllRSMs,3),RSMSize,'single');
    for rep = 1:RSMSize
        tempSubSample=SubSample;
        tempSubSample(rep,1)=0;
        tempRSMs=zeros(SubRSMVertSize,size(AllRSMs,2),size(AllRSMs,3),'single');
        for j=1:size(AllRSMs,3)
            [ ~,tempRSMs(:,:,j) ] = SubSampleVertRSM2SymRSM( AllRSMs(:,:,j),tempSubSample);
        end
        if ndims(AllRSMs) == 3
            [ OutVars ] = FullRCA( tempRSMs,'vertIn',1,...
                'DownSample',DownSample,...
                'DownSampleN',DownSampleN,...
                'DownSampleReps',DownSampleReps,...
                'RunPermute',RunPermute,...
                'PermuteReps',PermuteReps,...
                'OutputPermDists',OutputPermDists,...
                'zNorm',zNormRSM);
        else
            [ OutVars ] = FullRCA( tempRSMs,'vertIn',0,...
                'DownSample',DownSample,...
                'DownSampleN',DownSampleN,...
                'DownSampleReps',DownSampleReps,...
                'RunPermute',RunPermute,...
                'PermuteReps',PermuteReps,...
                'OutputPermDists',OutputPermDists,...
                'zNorm',zNormRSM);            
        end
%             if rep==1
%                 RFMats=single(OutVars.RF)/1000;
%                 RCMats=single(OutVars.RCMat)/1000;
%             else    
%                 RFMats=RFMats+single(OutVars.RF)/1000;
%                 RCMats=RCMats+single(OutVars.RCMat)/1000; 
%             end
        ALLRF(:,:,rep)=single(OutVars.RF);
        ALLRC(:,:,:,rep)=single(OutVars.RCMat);       
    end
    ALLRFDiff=normalize(ALLRF,3,'center')*-1;
    ALLRCDiff=normalize(ALLRC,4,'center')*-1;
    VertLabels=labels2uppertriuvectorlabels( UseLabels );
    for loadInd = 1:length(UseInd)
        loadNum=UseInd(loadInd,1);
        RCA_RF=squeeze(ALLRF(:,loadInd,:))';
        RCA_RF_Diff=squeeze(ALLRFDiff(:,loadInd,:))';
        try
            RCA_RF=array2table(RCA_RF,'VariableNames',UseLabels,'RowNames',ConditionNames);
            RCA_RF_Diff=array2table(RCA_RF_Diff,'VariableNames',UseLabels,'RowNames',ConditionNames);
        catch
            for i = 1:length(UseLabels)
                UseLabels{i,1}=[UseLabels{i,1},num2str(i)];
            end
            RCA_RF=array2table(RCA_RF,'VariableNames',UseLabels,'RowNames',ConditionNames);
            RCA_RF_Diff=array2table(RCA_RF_Diff,'VariableNames',UseLabels,'RowNames',ConditionNames);
        end            
        RCA_RC=mat2uppertriuvectormat(squeeze(ALLRC(:,:,loadInd,:)))';
        RCA_RC_Diff=mat2uppertriuvectormat(squeeze(ALLRCDiff(:,:,loadInd,:)))';           	
        try
            RCA_RC=array2table(RCA_RC,'VariableNames',VertLabels,'RowNames',ConditionNames);
            RCA_RC_Diff=array2table(RCA_RC_Diff,'VariableNames',VertLabels,'RowNames',ConditionNames);
        catch
            for i = 1:length(VertLabels)
                VertLabels{i,1}=[VertLabels{i,1},VertLabels(i)];
            end
            VertLabel=truncateStrCell(VertLabels);
            RCA_RC=array2table(RCA_RC,'VariableNames',VertLabels,'RowNames',ConditionNames);
            RCA_RC_Diff=array2table(RCA_RC_Diff,'VariableNames',VertLabels,'RowNames',ConditionNames);
        end     
        if ~exist(SaveDirs_RCA.(parcelName){loadNum,1},'file')
            mkdir(SaveDirs_RCA.(parcelName){loadNum,1});
        end
        SaveName=[SaveDirs_RCA.(parcelName){loadNum,1},SaveFileNames_RCA.(parcelName){loadNum,1}];
%         [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(RCA_RC),'ThresholdType','rank','Threshold',0.5);
%         rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
%         for nameNum=1:length(rowNames50)
%             rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
%         end
%         tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);
%         
%         [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(RCA_RC),'ThresholdType','rank','Threshold',0.25);
%         rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
%         for nameNum=1:length(rowNames25)
%             rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
%         end        
%         tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);
%         
%         [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(RCA_RC),'ThresholdType','rank','Threshold',0.10);
%         rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
%         for nameNum=1:length(rowNames10)
%             rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
%         end        
%         tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);
%         
%         [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(RCA_RC),'ThresholdType','rank','Threshold',0.05);
%         rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
%         for nameNum=1:length(rowNames05)
%             rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
%         end        
%         tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);
%         try
%             RCA_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
%             RCA_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
%             RCA_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
%             IndvRCA_RC=array2table(corrcoef(AllRSMs(:,:,loadInd)),'VariableNames',UseLabels,'RowNames',UseLabels);
%             RCA_RC_Z=array2table(vertRSM2SymRSM(nan_zscore(mat2uppertriuvectormat(table2array(RCA_RC))),[],1),'VariableNames',RCA_RC.Properties.VariableNames,'RowNames',RCA_RC.Properties.RowNames);
%             IndvRCA_RC_Z=array2table(vertRSM2SymRSM(nan_zscore(mat2uppertriuvectormat(table2array(IndvRCA_RC))),[],1),'VariableNames',IndvRCA_RC.Properties.VariableNames,'RowNames',IndvRCA_RC.Properties.RowNames);
%             RCA_RF_Z = array2table(nan_zscore(table2array(RCA_RF)')','VariableNames',RCA_RF.Properties.VariableNames);
%             skipGroupRC=0;
%         catch
%             RCA_RF_Z = array2table(nan_zscore(table2array(RCA_RF)')','VariableNames',RCA_RF.Properties.VariableNames);
%             RCA_RC_WholeCM=[];
%             RCA_RC_Degree=[];
%             IndvRCA_RC=[];
%             IndvRCA_RC_Z=[];
%             RCA_RC_Z=[];
%             skipGroupRC=1;
%         end
%             
%         if OutputPermDists == 1
%             RF_PermDist=OutVars.RF_Dist(:,:,loadInd);
%             RCMat_PermDist=mat2uppertriuvectormat(vertRSM2SymRSM(OutVars.RC_Dist(:,:,loadInd)'))';
%             save(SaveName,'RCA_RF','RCA_RC','IndvRCA_RC','RCA_RC_Degree','RF_PermDist','RCMat_PermDist','AnalysisParameters');
%         else
        save(SaveName,'RCA_RF','RCA_RC','RCA_RF_Diff','RCA_RC_Diff','AnalysisParameters');
     
    end
    %% compute parcellation and SS-based RSMs for 3rd order RSA.
    
%     ssRSM_RF=single(corrcoef(RFMats));
%     parcelRSM_RF=single(corrcoef(RFMats'));
%     ssRSM_RF_scaled=single(corrcoef(scaleVals(RFMats,2)));
%     parcelRSM_RF_scaled=single(corrcoef(scaleVals(RFMats',2))); 
%     if skipGroupRC==0
%         RCMatsVert=mat2uppertriuvectormat(RCMats);
%         if length(UseLabels) <= 200
%             ssRSM_RC=single(corrcoef(RCMatsVert));
%             parcelRSM_RC=single(corrcoef(RCMatsVert'));
%             ssRSM_RC_scaled=single(corrcoef(scaleVals(RCMatsVert,2)));
%             parcelRSM_RC_scaled=single(corrcoef(scaleVals(RCMatsVert',2)));  
%         else
%             ssRSM_RC=[];
%             parcelRSM_RC=[];
%             ssRSM_RC_scaled=[];
%             parcelRSM_RC_scaled=[];
%         end
% 
%         [Group_RF,Group_RC,Group_RCMat_Ts,Group_RCMat_Means,RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(RFMats,RCMats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
%             'vecName','RF','matName','RC','permZMat',OutVars.RCMat_permZ,'permPMat',OutVars.RCMat_permP,'permZVec',OutVars.RF_permZ,'permPVec',OutVars.RF_permP,'permZMat_Corrected',OutVars.RCMat_permZ_Corrected,'permPMat_Corrected',OutVars.RCMat_permP_Corrected,'permZVec_Corrected',OutVars.RF_permZ_Corrected,'permPVec_Corrected',OutVars.RF_permP_Corrected);
%         if ComputeReliability == 1
%             RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
%             RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
%             RF_losoReliability=ReliabilityResults.Vec_losoReliability;
%             RC_losoReliability=ReliabilityResults.Mat_losoReliability;
%             save([GroupAnalysisDir,'Group_ParcelReliability'],'RF_SplitHalfReliability','RC_SplitHalfReliability','RF_losoReliability','RC_losoReliability'); 
%         end
%         ssRSM_RSA = zeros(height(UseTable),height(UseTable),length(UseLabels),'single');
%         for j = 1:length(UseLabels)
%             ssRSM_RSA(:,:,j)=corrcoef(squeeze(AllRSMs(:,j,:)));
%         end   
%         ssRSM_RSA=array2table(mat2uppertriuvectormat(ssRSM_RSA),'VariableNames',UseLabels);
%         save([GroupAnalysisDir,'Group_RC'],'Group_RC');
%         save([GroupAnalysisDir,'Group_RCMat'],'Group_RCMat_Ts','Group_RCMat_Means','RC_Degree_Table');  
%         save([GroupAnalysisDir,'ThirdOrderRSMs_RC'],'ssRSM_RC','parcelRSM_RC','ssRSM_RC_scaled','parcelRSM_RC_scaled','ssRSM_RSA','UseInd','UseTable');
%     else
%         [Group_RF,ReliabilityResults] = GroupVecParcellationSummaryFigures(RFMats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'vecName','RF');
%         if ComputeReliability == 1
%             RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
%             RF_losoReliability=ReliabilityResults.Vec_losoReliability;
%             save([GroupAnalysisDir,'Group_ParcelReliability'],'RF_SplitHalfReliability','RF_losoReliability'); 
%         end   
%     end

    [GroupRCDiff_T,GroupRCDiff_Z,GroupRCDiff_M]=getTval(ALLRCDiff,3);
    GroupRCDiff_T=array2table(mat2uppertriuvectormat(squeeze(GroupRCDiff_T))','VariableNames',VertLabels,'RowNames',ConditionNames);
    GroupRCDiff_Z=array2table(mat2uppertriuvectormat(squeeze(GroupRCDiff_Z))','VariableNames',VertLabels,'RowNames',ConditionNames);
    GroupRCDiff_M=array2table(mat2uppertriuvectormat(squeeze(GroupRCDiff_M))','VariableNames',VertLabels,'RowNames',ConditionNames);
    GroupRC_M=array2table(mat2uppertriuvectormat(squeeze(nanmean(ALLRC,3)))','VariableNames',VertLabels,'RowNames',ConditionNames);

    [GroupRFDiff_T,GroupRFDiff_Z,GroupRFDiff_M]=getTval(ALLRFDiff,2);
    GroupRFDiff_T=array2table(squeeze(GroupRFDiff_T)','VariableNames',UseLabels,'RowNames',ConditionNames);
    GroupRFDiff_Z=array2table(squeeze(GroupRFDiff_Z)','VariableNames',UseLabels,'RowNames',ConditionNames);
    GroupRFDiff_M=array2table(squeeze(GroupRFDiff_M)','VariableNames',UseLabels,'RowNames',ConditionNames);
    GroupRF_M=array2table(squeeze(nanmean(ALLRF,2))','VariableNames',UseLabels,'RowNames',ConditionNames);

    save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters','UseInd','UseTable');
    save([GroupAnalysisDir,'Group_RF_ByItem'],'GroupRFDiff_T','GroupRFDiff_Z','GroupRFDiff_M','GroupRF_M');
    try
        save([GroupAnalysisDir,'Group_RC_ByItem'],'GroupRCDiff_T','GroupRCDiff_Z','GroupRCDiff_M','GroupRC_M');
    catch
        save([GroupAnalysisDir,'Group_RC_ByItem'],'GroupRCDiff_T','GroupRCDiff_Z','GroupRCDiff_M','GroupRC_M','-v7.3');
    end 
    toc
end
end 

