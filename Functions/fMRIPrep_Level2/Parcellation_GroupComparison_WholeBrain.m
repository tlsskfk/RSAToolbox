function [AnalysisParameters] = Parcellation_GroupComparison_WholeBrain(ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
%Overwrite previously saved files (default is no or 0; yes = 1)
if nargin==0    
    ExperimentsDir=uigetdir('/','Select Experiments Directory');
    ExperimentsDir=[ExperimentsDir,'/'];
end

[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Cell containing list of table names
[fmriprep_table_names] = VariableSetter('fmriprep_table_names',[],varargin);
%Cell containing fmriprep_tables to compare
[fmriprep_tables] = VariableSetter('fmriprep_tables',[],varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Subject or run level analysis. Will prompt request.
[MapAnalysisType] = VariableSetter('MapAnalysisType',[],varargin);
[MapAnalysisName] = VariableSetter('MapAnalysisName',[],varargin);
[MapVarName] = VariableSetter('MapVarName',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
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
    fmriprep_table_names=AnalysisParameters.fmriprep_table_names;
    MapVarName=AnalysisParameters.MapVarName;
    UncorrectedThresholds=AnalysisParameters.UncorrectedThresholds;
    ClusterThreshReps=AnalysisParameters.ClusterThreshReps;
else
    AnalysisParameters=struct;

end

if isempty(fmriprep_table_names)
    NumGroups=uiEnterName('2','Enter # of groups to compare');
    NumGroups=str2num(NumGroups);
    for i = 1:NumGroups
        [~,~,fmriprep_table_names{i,1}] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
    end
    fmriprep_table_names=fmriprep_table_names(:);
    [~,NewVarNames] = uiTableVarNames(fmriprep_table_names);
    fmriprep_table_names=[fmriprep_table_names,NewVarNames(:)];
end

fmriprep_tables=[];
GroupVar=[];
ComparisonNames=[];
NumGroups=size(fmriprep_table_names,1);
for i = 1:NumGroups
    [~,tempTable] = load_fmriprep_table('ExperimentsDir',ExperimentsDir,'fmriprep_table_name',fmriprep_table_names{i,1});
    if ~isempty(fmriprep_tables)
        tempVars=intersect(tempTable.Properties.VariableNames,fmriprep_tables.Properties.VariableNames,'stable');
        fmriprep_tables=fmriprep_tables(:,tempVars);
        tempTable=tempTable(:,tempVars); 
    end
    fmriprep_tables=[fmriprep_tables;tempTable];
    GroupVar=[GroupVar;ones(height(tempTable),1)*i];
    ComparisonNames=[ComparisonNames,fmriprep_table_names{i,2}];
end
GroupVar=array2table(GroupVar,'VariableNames',{'GroupVar'});
fmriprep_table=[fmriprep_tables,GroupVar];

AnalysisParameters.fmriprep_table_names=fmriprep_table_names;
AnalysisParameters.fmriprep_table=fmriprep_table;
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupDiff_GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
BrainMapDir=strrep(GroupDir,'/GroupDiff_GroupAnalysis/','/GroupDiff_GroupBrainMaps/');


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

if isempty(MapVarName)
   MapVarName=cell(NumGroups,1);
end

if isempty(MapAnalysisType)
   MapAnalysisType=cell(NumGroups,1);
end
if isempty(MapAnalysisName)
   MapAnalysisName=cell(NumGroups,1);
end
filePaths_MapData=cell(NumGroups,1);
SearchlightNames=cell(NumGroups,1);

%% Compile filepaths for input files for the analysis
for i = 1:NumGroups
    [filePaths_MapData{i,1},MapAnalysisType{i,1},MapAnalysisName{i,1}] = BIDsDirSearch(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),'SubjectOrRun',SubjectOrRun,'AnalysisType',MapAnalysisType{i,1},'AnalysisName',MapAnalysisName{i,1},'TitleTextName',['Select source of Map values for',newline,fmriprep_table_names{i,2}]);
    [SearchlightNames{i,1}] = uiNameSelect(filePaths_MapData{i,1}.Properties.VariableNames,['Select searchlight map for',newline,fmriprep_table_names{i,2}] ,1);
    if iscell(SearchlightNames{i,1})
        SearchlightNames{i,1}=SearchlightNames{i,1}{1,1};
    end
    if isempty(MapVarName{i,1})
        tempMapVarName=[];
        for j = 1:size(filePaths_MapData{i,1},1)
            try
                variableInfo = who('-file', filePaths_MapData{i,1}.(SearchlightNames{i,1}){j,1});
                tempMapVarName=unique([tempMapVarName;variableInfo]);
            catch
                continue
            end
        end
        [MapVarName{i,1}] = uiNameSelect(tempMapVarName,['Select Map variable for',newline,fmriprep_table_names{i,2}] ,1);    
    end
end
AnalysisParameters.MapAnalysisType=MapAnalysisType;
AnalysisParameters.MapAnalysisName=MapAnalysisName;
AnalysisParameters.MapVarName=MapVarName;   

AnalysisType=['GroupDiff_',MapAnalysisType{1,1}];
AnalysisNamesSuffix=[MapAnalysisName{1,1},'_',MapVarName{1,1}];
if DefaultName==0
    if isempty(AnalysisName)
        AnalysisName=uiEnterName(['GroupDiff_',AnalysisNamesSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    end
    ComparisonNames=uiEnterName([ComparisonNames],['Enter name for ',AnalysisType,newline,'analysis below:']);
else
    AnalysisName=['GroupDiff_',AnalysisNamesSuffix];
end
if bySS == 1
    GroupAnalysisDir=[GroupDir,ComparisonNames,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightNames{1,1},'/'];
    GroupBrainMapsDir=[BrainMapDir,ComparisonNames,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',SearchlightNames{1,1},'/'];
else
    GroupAnalysisDir=[GroupDir,ComparisonNames,'/',AnalysisType,'/',AnalysisName,'/',SearchlightNames{1,1},'/'];
    GroupBrainMapsDir=[BrainMapDir,ComparisonNames,'/',AnalysisType,'/',AnalysisName,'/',SearchlightNames{1,1},'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);     
end
if ~exist(GroupBrainMapsDir,'file')
    mkdir(GroupBrainMapsDir);     
end    


Use_fmriprep_table=[];
MapData=[];
for i = 1:NumGroups
    [tempMapDataAll,compiled_fmriprep_table,tempMapataSampleInd] = CompileND(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),...
        'LoadVarName',MapVarName{i,1},...
        'LoadVarFormat','Array',...
        'DataType','ByParcellation',...
        'AnalysisType',MapAnalysisType{i,1},...
        'AnalysisName',MapAnalysisName{i,1},...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Array',...
        'ParcelName',SearchlightNames{i,1});  
    Use_fmriprep_table=[Use_fmriprep_table;compiled_fmriprep_table];
    MapData=cat(4,MapData,tempMapDataAll);     
end
GroupVar=Use_fmriprep_table.GroupVar;

DiffMap_Mean=nanmean(MapData(:,:,:,GroupVar==1),4)-nanmean(MapData(:,:,:,GroupVar==2),4);
[~,pMap]=ttest2(MapData(:,:,:,GroupVar==1),MapData(:,:,:,GroupVar==2),'dim',4);
DiffMap_Z=p2z(pMap,2);
DiffMap_Z(isinf(DiffMap_Z))=nan;
DiffMap_Z=DiffMap_Z.*sign(DiffMap_Mean);

SaveBrik_3mmMNI(cat(4,DiffMap_Mean,DiffMap_Z),{'DiffMap_Mean','DiffMap_Z'},[GroupBrainMapsDir,'DiffMap_SearchlightMaps']);
save([GroupAnalysisDir,'DiffMap_SearchlightMaps'],'DiffMap_Mean','DiffMap_Z','AnalysisParameters','Use_fmriprep_table');

GroupClusterDists=nan(ClusterThreshReps,length(UncorrectedThresholds));
for rep =  1:ClusterThreshReps
    disp(['Running cluster correction: ',num2str((rep/ClusterThreshReps)*100),'% Complete.'])    
    GroupVar=GroupVar(randperm(length(GroupVar)));
    tempDiffMap_Mean=nanmean(MapData(:,:,:,GroupVar==1),4)-nanmean(MapData(:,:,:,GroupVar==2),4);
    [~,pMap]=ttest2(MapData(:,:,:,GroupVar==1),MapData(:,:,:,GroupVar==2),'dim',4);
    tempDiffMap_Z=p2z(pMap,2);
    tempDiffMap_Z(isinf(tempDiffMap_Z))=nan;
    tempDiffMap_Z=tempDiffMap_Z.*sign(tempDiffMap_Mean);
    [maskcoords,vals]=mat2coords(tempDiffMap_Z);        
    [ Clusters ] = ClusterFinder(vals,maskcoords,'both','stat',UncorrectedThresholds,1,'side');
    GroupClusterDists(rep,:)=Clusters.MaxClusterSizes;    
end

GroupClusterThreshs=prctile(GroupClusterDists,95);    
[maskcoords,vals]=mat2coords(DiffMap_Z);  
ThreshDiffMap_Z=[];
mapLabels=[];
Clusters=cell(length(UncorrectedThresholds),1);
for j = 1:length(UncorrectedThresholds)
    [ Clusters{j,1} ] = ClusterFinder(vals,maskcoords,'both','stat',UncorrectedThresholds(1,j),GroupClusterThreshs(1,j),'side');
    if ~isempty( Clusters{j,1}.ThresholdedClusters{1,1}{1,1})
        CorrectedCoords=[];
        for k = 1:size(Clusters{j,1}.ThresholdedClusters{1,1},1)
            CorrectedCoords=[CorrectedCoords;Clusters{j,1}.ThresholdedClusters{1,1}{k,1}];
        end
        threshmat=coords2mat(CorrectedCoords,DiffMap_Z,ones(size(CorrectedCoords,1),1));
        mapLabels=[mapLabels,{['VoxThresh_Zgt',num2str4filename(UncorrectedThresholds(1,j),2)]}];
        ThreshDiffMap_Z=cat(4,ThreshDiffMap_Z,DiffMap_Z.*threshmat);
    end
end
if ~isempty(ThreshDiffMap_Z)
    SaveBrik_3mmMNI(ThreshDiffMap_Z,mapLabels,[GroupBrainMapsDir,'DiffMap_Z_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)]);
    save([GroupAnalysisDir,'DiffMap_Z_CorrectedSLMaps_PermReps',num2str(ClusterThreshReps)],'ThreshDiffMap_Z','GroupClusterDists','GroupClusterThreshs','UncorrectedThresholds','Clusters','AnalysisParameters');     
end
end
