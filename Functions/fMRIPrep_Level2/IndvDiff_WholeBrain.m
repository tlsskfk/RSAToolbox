function [AnalysisParameters] = IndvDiff_WholeBrain(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
%Overwrite previously saved files (default is no or 0; yes = 1)
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end

[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Subject or run level analysis. Will prompt request.
[MapAnalysisType] = VariableSetter('MapAnalysisType',[],varargin);
[MapAnalysisName] = VariableSetter('MapAnalysisName',[],varargin);
[corrType] = VariableSetter('corrType',[],varargin);
[MapVarName] = VariableSetter('MapVarName',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[IncludeTaskVars] = VariableSetter('IncludeTaskVars',[],varargin); %1 = yes 
[TaskVarNames] = VariableSetter('TaskVarNames',[],varargin);
[IncludeSSVars] = VariableSetter('IncludeSSVars',[],varargin);
[SSVarNames] = VariableSetter('SSVarNames',[],varargin);
[zNorm] = VariableSetter('zNorm',[],varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
[UncorrectedThresholds] = VariableSetter('UncorrectedThresholds',[1.96,2.576,2.807,3.2905],varargin);
[ClusterThreshReps] = VariableSetter('ClusterThreshReps',1000,varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            ExperimentsDir=strrep(ExperimentsDir,'\','/');
            ParamDir=uigetdir(ExperimentsDir);
            ParamDir=[ParamDir,'/'];
            [~,TableTypeFileNames]=getFolderAndFileNames(ParamDir);
            TableTypeFileNames=uiNameSelect(TableTypeFileNames,'Select table file',1); 
            if iscell(TableTypeFileNames)
                TableTypeFileNames=TableTypeFileNames{1,1};
            end   
            TableTypeAnalysisParams=load([ParamDir,TableTypeFileNames],'AnalysisParameters');
            AnalysisParameters=TableTypeAnalysisParams.AnalysisParameters;
        catch
            disp('No Existing Parameters!')
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
    MapAnalysisType=AnalysisParameters.MapAnalysisType;
    MapAnalysisName=AnalysisParameters.MapAnalysisName;
    MapVarName=AnalysisParameters.MapVarName;
    UncorrectedThresholds=AnalysisParameters.UncorrectedThresholds;
    ClusterThreshReps=AnalysisParameters.ClusterThreshReps;
else
    AnalysisParameters=struct;

end

ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','IndvDiff_GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
BrainMapDir=strrep(GroupDir,'/IndvDiff_GroupAnalysis/','/IndvDiff_GroupBrainMaps/');


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
AnalysisParameters.SubjectOrRun=SubjectOrRun;
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=dataInd;
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end

%% Compile filepaths for input files for the analysis

[filePaths_MapData,MapAnalysisType,MapAnalysisName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType',MapAnalysisType,'AnalysisName',MapAnalysisName,'TitleTextName',['Select source of Map values']);
[SearchlightNames] = uiNameSelect(filePaths_MapData.Properties.VariableNames,['Select searchlight map'] ,1);
if iscell(SearchlightNames)
    SearchlightNames=SearchlightNames{1,1};
end
if iscell(MapAnalysisType)
    MapAnalysisType=MapAnalysisType{1,1};
end
if iscell(MapAnalysisName)
    MapAnalysisName=MapAnalysisName{1,1};
end

if isempty(MapVarName)
    tempMapVarName=[];
    for j = 1:size(filePaths_MapData,1)
        try
            variableInfo = who('-file', filePaths_MapData.(SearchlightNames){j,1});
            tempMapVarName=unique([tempMapVarName;variableInfo]);
        catch
            continue
        end
    end
    [MapVarName] = uiNameSelect(tempMapVarName,['Select Map variable'] ,1);    
end
zSuffix=[];
if isempty(zNorm)
    SingleSelect=1; %Allows only a single value to be selected.
    [zNorm] = uiNameSelect({'Yes','No'},'zNorm data?',SingleSelect); 
    if strcmp(zNorm,'Yes')
        zNorm=1;
    else
        zNorm=0;
    end
end
if zNorm==1
    zSuffix='_zNormed';
end
if isempty(corrType)
    SingleSelect=1; %Allows only a single value to be selected.
    [corrType] = uiNameSelect({'Pearson','Spearman','Kendall'},'Select correlation type:',SingleSelect);
end
AnalysisParameters.MapAnalysisType=MapAnalysisType;
AnalysisParameters.MapAnalysisName=MapAnalysisName;
AnalysisParameters.MapVarName=MapVarName;   
AnalysisParameters.zNorm=zNorm; 
AnalysisType=['IndvDiff_',MapAnalysisType];
AnalysisNamesSuffix=[MapAnalysisName,'_',MapVarName,zSuffix,'_',corrType];
if DefaultName==0
    if isempty(AnalysisName)
        AnalysisName=uiEnterName(['IndvDiff_',AnalysisNamesSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    end
else
    AnalysisName=['IndvDiff_',AnalysisNamesSuffix];
end
if bySS == 1
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightNames,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightNames,'/'];
else
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SearchlightNames,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',SearchlightNames,'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);     
end
if ~exist(GroupBrainMapsDir,'file')
    mkdir(GroupBrainMapsDir);     
end    

[BrainData,compiled_fmriprep_table] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName',MapVarName,...
        'LoadVarFormat','Array',...
        'DataType','ByParcellation',...
        'AnalysisType',MapAnalysisType,...
        'AnalysisName',MapAnalysisName,...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Array',...
        'ParcelName',SearchlightNames);  
BrainSize=size(BrainData);
if length(BrainSize)==5
    [ConditionNames] = CompileND(ExperimentsDir,fmriprep_table,...
            'LoadVarName','ConditionNames',...
            'LoadVarFormat','Array',...
            'DataType','ByParcellation',...
            'AnalysisType',MapAnalysisType,...
            'AnalysisName',MapAnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'ParcelName',SearchlightNames); 
    ConditionNames=ConditionNames(:,1)';
    ConditionName=uiNameSelect(ConditionNames,'Select glm condition',1);
    BrainData=squeeze(BrainData(:,:,:,ismember(ConditionNames,ConditionName),:));

else
    ConditionNames='';
    ConditionName='';
end
BrainSize=size(BrainData);
BrainData=reshape(BrainData,[BrainSize(1,1)*BrainSize(1,2)*BrainSize(1,3),BrainSize(1,4)])';
selectInds=~isnan(BrainData);
groupingInds=sum(single(selectInds),1);
useInds=groupingInds~=0;
groupingInds=groupingInds(:,useInds);
fullSetInd=groupingInds==size(BrainData,1);
BrainData=BrainData(:,useInds);
[UseCoords]=mat2coords(reshape(useInds,[BrainSize(1,1),BrainSize(1,2),BrainSize(1,3)]));
numVoxels=size(BrainData,2);
if zNorm==1
    BrainData=nan_zscore(BrainData')';
end
SingleSelect=1; %Allows only a single value to be selected.
AllVarNames=[];
if isempty(IncludeTaskVars)
    [IncludeTaskVars] = uiNameSelect({'Yes','No'},'Include task vars?',SingleSelect);
    if strcmpi(IncludeTaskVars,'Yes')
        IncludeTaskVars=1;
    else
        IncludeTaskVars=0;
    end
end
if IncludeTaskVars==1
    [TaskVars_data_all,~,TaskVarInd] = CompileND(ExperimentsDir,compiled_fmriprep_table,...
        'LoadVarName','beh_vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','beh_vars',...
        'AnalysisName','Overall',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    if isempty(TaskVarNames)
        TaskVarNames=uiNameSelect(TaskVars_data_all.Properties.VariableNames,'Select task varaibles to compile.',0);
    end
    TaskVars_data_all=TaskVars_data_all(:,TaskVarNames);
    AllVarNames=TaskVarNames(:);
end
%%Compile SS variables
SingleSelect=1; %Allows only a single value to be selected.
if isempty(IncludeSSVars)
    [IncludeSSVars] = uiNameSelect({'Yes','No'},'Include SS vars?',SingleSelect);
    if strcmpi(IncludeSSVars,'Yes')
        IncludeSSVars=1;
    else
        IncludeSSVars=0;
    end
end
if IncludeSSVars==1
    [SSVars_data_all,~,SSVarInd] = CompileND(ExperimentsDir,compiled_fmriprep_table,...
        'LoadVarName','SS_Vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','SSVars',...
        'AnalysisName','SSVars_All',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    if isempty(SSVarNames)
        SSVarNames=uiNameSelect(SSVars_data_all.Properties.VariableNames,'Select SS varaibles to compile.',0);
    end
    SSVars_data_all=SSVars_data_all(:,SSVarNames);
    AllVarNames=[AllVarNames;SSVarNames(:)];
end
ClusterThreshInfo=cell(size(AllVarNames,1),3);
ClusterThreshInfo(:,1)=AllVarNames;
BrainMaps=[];
CorrectedBrainMaps=[];
BrainLabels=[];
CorrectedBrainLabels=[];
count=1;
if IncludeTaskVars==1
    UseBrainData=BrainData(TaskVarInd==1,:);
    for i = 1:length(TaskVarNames)
        tic
        GroupClusterDists=nan(ClusterThreshReps,length(UncorrectedThresholds));
        AllRs=nan(1,numVoxels);
        AllPs=nan(1,numVoxels);           
        TaskVarName=TaskVarNames{i};
        disp(TaskVarName);
        behData=TaskVars_data_all.(TaskVarName);
        nanInd=~isnan(behData);
        behData=behData(nanInd,:);
        tempBrainData=UseBrainData(nanInd,:);
        [AllRs(1,fullSetInd==1),AllPs(1,fullSetInd==1)]=corr(behData,tempBrainData(:,fullSetInd==1),'type',corrType);
        [AllRs(1,fullSetInd==0),AllPs(1,fullSetInd==0)]=corr(behData,tempBrainData(:,fullSetInd==0),'rows','pairwise','type',corrType);
        AllZs=p2z(AllPs,2).*sign(AllRs);
        zMap=coords2mat(UseCoords,zeros(BrainSize(1,1),BrainSize(1,2),BrainSize(1,3)),AllZs(:));
        rMap=coords2mat(UseCoords,zeros(BrainSize(1,1),BrainSize(1,2),BrainSize(1,3)),AllRs(:));
        BrainMaps=cat(4,BrainMaps,rMap);
        BrainMaps=cat(4,BrainMaps,zMap);
        BrainLabels=[BrainLabels,{[TaskVarName,'_r'],[TaskVarName,'_Z']}];
        parfor j = 1:ClusterThreshReps
            tempPs=AllPs*0;
            tempRs=AllRs*0;
            tempbehData=behData(randperm(length(behData)));
            [tempRs(1,fullSetInd==1),tempPs(1,fullSetInd==1)]=corr(tempbehData,tempBrainData(:,fullSetInd==1),'type',corrType);
            [tempRs(1,fullSetInd==0),tempPs(1,fullSetInd==0)]=corr(tempbehData,tempBrainData(:,fullSetInd==0),'rows','pairwise','type',corrType);                
            tempZs=p2z(tempPs,2).*sign(tempRs);
            [ Clusters ] = ClusterFinder(tempZs',UseCoords,'both','stat',UncorrectedThresholds,1,'side');
            GroupClusterDists(j,:)=Clusters.MaxClusterSizes;                    
        end
        GroupClusterThreshs=prctile(GroupClusterDists,95); 
        ClusterThreshInfo{count,2}=GroupClusterThreshs;
        ClusterThreshInfo{count,3}=GroupClusterDists;
        Clusters=cell(length(UncorrectedThresholds),1);
        for j = 1:length(UncorrectedThresholds)
            [ Clusters{j,1} ] = ClusterFinder(AllZs',UseCoords,'both','stat',UncorrectedThresholds(1,j),GroupClusterThreshs(1,j),'side');
            if ~isempty( Clusters{j,1}.ThresholdedClusters{1,1}{1,1})
                CorrectedCoords=[];
                for k = 1:size(Clusters{j,1}.ThresholdedClusters{1,1},1)
                    CorrectedCoords=[CorrectedCoords;Clusters{j,1}.ThresholdedClusters{1,1}{k,1}];
                end
                threshmat=coords2mat(CorrectedCoords,zMap,ones(size(CorrectedCoords,1),1));
                CorrectedBrainLabels=[CorrectedBrainLabels,{[TaskVarName,'_Zgt',num2str4filename(UncorrectedThresholds(1,j),2)]}];
                CorrectedBrainMaps=cat(4,CorrectedBrainMaps,zMap.*threshmat);
            end
        end
        count=count+1;
        toc
    end
end

if IncludeSSVars==1
    UseBrainData=BrainData(SSVarInd==1,:);
    for i = 1:length(SSVarNames)
        tic
        GroupClusterDists=nan(ClusterThreshReps,length(UncorrectedThresholds));
        AllRs=nan(1,numVoxels);
        AllPs=nan(1,numVoxels);           
        SSVarName=SSVarNames{i};
        disp(SSVarName);
        behData=SSVars_data_all.(SSVarName);
        nanInd=~isnan(behData);
        behData=behData(nanInd,:);
        tempBrainData=UseBrainData(nanInd,:);
        [AllRs(1,fullSetInd==1),AllPs(1,fullSetInd==1)]=corr(behData,tempBrainData(:,fullSetInd==1),'type',corrType);
        [AllRs(1,fullSetInd==0),AllPs(1,fullSetInd==0)]=corr(behData,tempBrainData(:,fullSetInd==0),'rows','pairwise','type',corrType);
        AllZs=p2z(AllPs,2).*sign(AllRs);
        zMap=coords2mat(UseCoords,zeros(BrainSize(1,1),BrainSize(1,2),BrainSize(1,3)),AllZs(:));
        rMap=coords2mat(UseCoords,zeros(BrainSize(1,1),BrainSize(1,2),BrainSize(1,3)),AllRs(:));
        BrainMaps=cat(4,BrainMaps,rMap);
        BrainMaps=cat(4,BrainMaps,zMap);
        BrainLabels=[BrainLabels,{[SSVarName,'_r'],[SSVarName,'_Z']}];
        parfor j = 1:ClusterThreshReps
            tempPs=AllPs*0;
            tempRs=AllRs*0;
            tempbehData=behData(randperm(length(behData)));
            [tempRs(1,fullSetInd==1),tempPs(1,fullSetInd==1)]=corr(tempbehData,tempBrainData(:,fullSetInd==1),'type',corrType);
            [tempRs(1,fullSetInd==0),tempPs(1,fullSetInd==0)]=corr(tempbehData,tempBrainData(:,fullSetInd==0),'rows','pairwise','type',corrType);                
            tempZs=p2z(tempPs,2).*sign(tempRs);
            [ Clusters ] = ClusterFinder(tempZs',UseCoords,'both','stat',UncorrectedThresholds,1,'side');
            GroupClusterDists(j,:)=Clusters.MaxClusterSizes;                    
        end
        GroupClusterThreshs=prctile(GroupClusterDists,95);
        ClusterThreshInfo{count,2}=GroupClusterThreshs;
        ClusterThreshInfo{count,3}=GroupClusterDists;        
        Clusters=cell(length(UncorrectedThresholds),1);
        for j = 1:length(UncorrectedThresholds)
            [ Clusters{j,1} ] = ClusterFinder(AllZs',UseCoords,'both','stat',UncorrectedThresholds(1,j),GroupClusterThreshs(1,j),'side');
            if ~isempty( Clusters{j,1}.ThresholdedClusters{1,1}{1,1})
                CorrectedCoords=[];
                for k = 1:size(Clusters{j,1}.ThresholdedClusters{1,1},1)
                    CorrectedCoords=[CorrectedCoords;Clusters{j,1}.ThresholdedClusters{1,1}{k,1}];
                end
                threshmat=coords2mat(CorrectedCoords,zMap,ones(size(CorrectedCoords,1),1));
                CorrectedBrainLabels=[CorrectedBrainLabels,{[SSVarName,'_Zgt',num2str4filename(UncorrectedThresholds(1,j),2)]}];
                CorrectedBrainMaps=cat(4,CorrectedBrainMaps,zMap.*threshmat);
            end
        end
        toc
        count=count+1;
    end
end

SaveBrik_3mmMNI(BrainMaps,BrainLabels,[GroupBrainMapsDir,'IndvDiff_SearchlightMaps',ConditionName]);
save([GroupAnalysisDir,'IndvDiff_SearchlightMaps',ConditionName],'BrainMaps','BrainLabels','AnalysisParameters','compiled_fmriprep_table');
save([GroupAnalysisDir,'IndvDiff_CorrectedSLMaps',ConditionName],'CorrectedBrainMaps','CorrectedBrainLabels','ClusterThreshInfo','UncorrectedThresholds','AnalysisParameters');     

if ~isempty(CorrectedBrainMaps)
    SaveBrik_3mmMNI(CorrectedBrainMaps,CorrectedBrainLabels,[GroupBrainMapsDir,'IndvDiff_CorrectedSLMaps',ConditionName]);
end

end
