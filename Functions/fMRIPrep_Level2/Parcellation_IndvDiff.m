function [AnalysisParameters] = Parcellation_IndvDiff(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
%Overwrite previously saved files (default is no or 0; yes = 1)
%Parcellation_IndvDiff(fmriprep_table,ExperimentsDir,'IncludeSSVars',1,'SSVarNames',,'beh_var_names',)
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
MakeFigs = VariableSetter('MakeFigs',[],varargin);
[AverageByVar] = VariableSetter('AverageByVar',[],varargin);
[AveragingVar] = VariableSetter('AveragingVar',[],varargin);
[SSVarNames] = VariableSetter('SSVarNames',[],varargin);
[IncludeSSVars] = VariableSetter('IncludeSSVars',[],varargin);
[beh_var_names] = VariableSetter('beh_var_names',[],varargin);
[motion_var_names] = VariableSetter('motion_var_names',[],varargin);
[NodeAnalysisType] = VariableSetter('NodeAnalysisType',[],varargin);
[NodeAnalysisName] = VariableSetter('NodeAnalysisName',[],varargin);
[EdgeAnalysisType] = VariableSetter('EdgeAnalysisType',[],varargin);
[EdgeAnalysisName] = VariableSetter('EdgeAnalysisName',[],varargin);
[RFX_Slopes] = VariableSetter('RFX_Slopes',[],varargin);
[RFX_Intercepts] = VariableSetter('RFX_Intercepts',[],varargin);
[RunVarTypes] = VariableSetter('RunVarTypes',[],varargin);
[NodeVarName] = VariableSetter('NodeVarName',[],varargin);
[EdgeVarName] = VariableSetter('EdgeVarName',[],varargin);
[NodeRowCompileVar] = VariableSetter('NodeRowCompileVar',[],varargin);
[EdgeRowCompileVar] = VariableSetter('EdgeRowCompileVar',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
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
    AverageByVar=AnalysisParameters.AverageByVar;
    AveragingVar=AnalysisParameters.AveragingVar;
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    MakeFigs=AnalysisParameters.MakeFigs;
    SSVarNames=AnalysisParameters.SSVarNames;
    IncludeSSVars=AnalysisParameters.IncludeSSVars;    
    NodeAnalysisType=AnalysisParameters.NodeAnalysisType;
    NodeAnalysisName=AnalysisParameters.NodeAnalysisName;
    EdgeAnalysisType=AnalysisParameters.EdgeAnalysisType;
    EdgeAnalysisName=AnalysisParameters.EdgeAnalysisName;
    beh_var_names=AnalysisParameters.beh_var_names;
    motion_var_names=AnalysisParameters.motion_var_names;
    fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    ParcelNames=AnalysisParameters.ParcelNames;   
    RFX_Slopes=AnalysisParameters.RFX_Slopes;
    RFX_Intercepts=AnalysisParameters.RFX_Intercepts;    
    RunVarTypes=AnalysisParameters.RunVarTypes;
    NodeVarName=AnalysisParameters.NodeVarName;
    EdgeVarName=AnalysisParameters.EdgeVarName;  
    try
        NodeRowCompileVar=AnalysisParameters.NodeRowCompileVar;
        EdgeRowCompileVar=AnalysisParameters.EdgeRowCompileVar;
    catch
        NodeRowCompileVar=[];
        EdgeRowCompileVar=[];
    end
else
    AnalysisParameters=struct;
    ParcelNames=[];
end

if isempty(fmriprep_table_name)
[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
AnalysisParameters.fmriprep_table_name=fmriprep_table_name;
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','IndvDiff_GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/IndvDiff_GroupAnalysis/','/IndvDiff_GroupFigures/');
BrainMapDir=strrep(GroupDir,'/IndvDiff_GroupAnalysis/','/IndvDiff_GroupBrainMaps/');

if isempty(fmriprep_table_name)
    [~,fmriprep_table_names]=getFolderAndFileNames('fmriprep_table/');
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect(fmriprep_table_names,'Select fmriprep_table used:',SingleSelect);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
end

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(MakeFigs)
    SingleSelect=1; %Allows only a single value to be selected.
    [MakeFigs] = uiNameSelect({'Yes','No'},'Make figures?',SingleSelect);
end
if strcmpi(MakeFigs,'Yes')
    MakeFigs=1;
else
    MakeFigs=0;
end
AnalysisParameters.MakeFigs=MakeFigs;
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
AveSuffix=[];
if isempty(AverageByVar)
    SingleSelect=1; %Allows only a single value to be selected.
    [AverageByVar] = uiNameSelect({'Yes','No'},'Collapse/average data prior to analysis?:',SingleSelect);    
    if strcmpi(AverageByVar,'Yes')
        SingleSelect=1;
        AveragingVar=uiNameSelect(fmriprep_table.Properties.VariableNames,'Select varaible to collapse/average over:',0);
        if ~iscell(AveragingVar)
            AveragingVar={AveragingVar};
        end
        AveSuffix=['_AveBy',AveragingVar{1,1}];
    else
        AveragingVar=[];
    end    
end

AnalysisParameters.AverageByVar=AverageByVar;
AnalysisParameters.AveragingVar=AveragingVar;
if isempty(RunVarTypes)
    SingleSelect=1; %Allows only a single value to be selected.
    [RunVarTypes] = uiNameSelect({'Nodes','Edges','Both'},'Run Indv Diff measures on parcel...',SingleSelect);
end
AnalysisParameters.RunVarTypes=RunVarTypes;


%% Compile filepaths for input files for the analysis
if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
    [filePaths_NodeData,NodeAnalysisType,NodeAnalysisName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType',NodeAnalysisType,'AnalysisName',NodeAnalysisName,'TitleTextName',['Select source of node values for',newline,'individual difference measures:']);
    AllParcelNamesNodes=filePaths_NodeData.Properties.VariableNames;
    if isempty(NodeVarName)
        for i = 1:size(filePaths_NodeData,1)
            for j = 1:length(AllParcelNamesNodes)
                try
                    variableInfo = who('-file', filePaths_NodeData.(AllParcelNamesNodes{1,j}){i,1});
                    NodeVarName=unique([NodeVarName;variableInfo]);
                catch
                    continue
                end
            end
        end
        [NodeVarName] = uiNameSelect(NodeVarName,'Select node variable to load:' ,1);    
    end
    if isempty(NodeRowCompileVar)
        for i = 1:size(filePaths_NodeData,1)
            for j = 1:length(AllParcelNamesNodes)
                try
                    tempNodeData=load(filePaths_NodeData.(AllParcelNamesNodes{1,j}){i,1},NodeVarName);
                    tempNodeData=tempNodeData.(NodeVarName);
                    if istable(tempNodeData)
                        if ~isempty(tempNodeData.Properties.RowNames)
                            NodeRowCompileVar=uiNameSelect(tempNodeData.Properties.RowNames,'Select node row variable:',1);
                            break
                        end
                    end    
                catch
                    continue
                end
            end
            if ~isempty(NodeRowCompileVar)
                break
            end
        end  
    elseif strcmpi(NodeRowCompileVar,'NA_EMPTY') 
        NodeRowCompileVar=[];
    end
else
    AllParcelNamesNodes=[];
    NodeAnalysisType=[];
    NodeAnalysisName=[];
    NodeVarName=[];
    NodeRowCompileVar=[];
end
AnalysisParameters.NodeAnalysisType=NodeAnalysisType;
AnalysisParameters.NodeAnalysisName=NodeAnalysisName;
AnalysisParameters.NodeVarName=NodeVarName;   
AnalysisParameters.NodeRowCompileVar=NodeRowCompileVar;   

if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
    [filePaths_EdgeData,EdgeAnalysisType,EdgeAnalysisName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType',EdgeAnalysisType,'AnalysisName',EdgeAnalysisName,'TitleTextName',['Select source of edge values for',newline,'individual difference measures:']);
    AllParcelNamesEdges=filePaths_EdgeData.Properties.VariableNames;
    if isempty(EdgeVarName)
        for i = 1:size(filePaths_EdgeData,1)
            for j = 1:length(AllParcelNamesEdges)
                try
                    variableInfo = who('-file', filePaths_EdgeData.(AllParcelNamesEdges{1,j}){i,1});
                    EdgeVarName=unique([EdgeVarName;variableInfo]);
                catch
                    continue
                end
            end
        end
        [EdgeVarName] = uiNameSelect(EdgeVarName,'Select edge variable to load:' ,1);
    end
    if isempty(EdgeRowCompileVar)
        for i = 1:size(filePaths_EdgeData,1)
            for j = 1:length(AllParcelNamesEdges)
                try
                    tempEdgeData=load(filePaths_EdgeData.(AllParcelNamesEdges{1,j}){i,1},EdgeVarName);
                    tempEdgeData=tempEdgeData.(EdgeVarName);
                    if istable(tempEdgeData)                    
                        if ~isempty(tempEdgeData.Properties.RowNames)
                            if ~any(~ismember(tempEdgeData.Properties.RowNames,tempEdgeData.Properties.VariableNames)) && length(tempEdgeData.Properties.RowNames)==length(tempEdgeData.Properties.VariableNames)
                                break
                            end
                            EdgeRowCompileVar=uiNameSelect(tempEdgeData.Properties.RowNames,'Select edge row variable:',1);
                            break
                        end
                    end    
                catch
                    continue
                end
            end
            if ~isempty(EdgeRowCompileVar)
                break
            end            
        end  
    elseif strcmpi(EdgeRowCompileVar,'NA_EMPTY') 
        EdgeRowCompileVar=[];
    end
else
    AllParcelNamesEdges=[];
    EdgeAnalysisType=[];
    EdgeAnalysisName=[];
    EdgeVarName=[];
    EdgeRowCompileVar=[];
end
AnalysisParameters.EdgeAnalysisType=EdgeAnalysisType;
AnalysisParameters.EdgeAnalysisName=EdgeAnalysisName;
AnalysisParameters.EdgeVarName=EdgeVarName; 
AnalysisParameters.EdgeRowCompileVar=EdgeRowCompileVar; 

if strcmpi(RunVarTypes,'Nodes')
    AllParcelNames=AllParcelNamesNodes;
    AnalysisType=['IndvDiff_',NodeAnalysisType];
    AnalysisNamesSuffix=[NodeAnalysisName,'_',NodeVarName,NodeRowCompileVar];
elseif strcmpi(RunVarTypes,'Edges')
    AllParcelNames=AllParcelNamesEdges;
    AnalysisType=['IndvDiff_',EdgeAnalysisType];
    AnalysisNamesSuffix=[EdgeAnalysisName,'_',EdgeVarName,EdgeRowCompileVar];   
elseif strcmpi(RunVarTypes,'Both')
    AllParcelNames=AllParcelNamesNodes(ismember(AllParcelNamesNodes,AllParcelNamesEdges)==1);
    AnalysisType=['IndvDiff_',NodeAnalysisType,'_',EdgeAnalysisType];
    if strcmpi(NodeAnalysisName,EdgeAnalysisName)
        AnalysisNamesSuffix=[NodeAnalysisName,'_',NodeVarName,NodeRowCompileVar,'_',EdgeVarName,EdgeRowCompileVar];
    else
        AnalysisNamesSuffix=[NodeAnalysisName,'_',EdgeAnalysisName,'_',NodeVarName,NodeRowCompileVar,'_',EdgeVarName,EdgeRowCompileVar];
    end
end
if DefaultName==0
    if isempty(AnalysisName)
        AnalysisName=uiEnterName(['IndvDiff_',AnalysisNamesSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    end
end
%%Compile Behavioral Data
try
[Behavior_data_all,Use_fmriprep_table,~] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','beh_vars',...
    'LoadVarFormat','Table',...
    'DataType','Other',...
    'AnalysisType','beh_vars',...
    'AnalysisName','Overall',...
    'SubjectOrRun',SubjectOrRun,...
    'TableOrArray','Table',...
    'ParcelName','None');
catch
[Behavior_data_all,Use_fmriprep_table,~] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','SS_Vars',...
    'LoadVarFormat','Table',...
    'DataType','Other',...
    'AnalysisType','SSVars',...
    'AnalysisName','SSVars_All',...
    'SubjectOrRun',SubjectOrRun,...
    'TableOrArray','Table',...
    'ParcelName','None');
    IncludeSSVars=0;
end
All_beh_var_names=Behavior_data_all.Properties.VariableNames;
if isempty(beh_var_names)
    beh_var_names=uiNameSelect(All_beh_var_names,'Select behavioral varaibles to compare.',0);
end
if ~iscell(beh_var_names)
    beh_var_names={beh_var_names};
end
AnalysisParameters.beh_var_names=beh_var_names;
Behavior_data_all=Behavior_data_all(:,beh_var_names);

if isempty(SSVarNames)
    SingleSelect=1; %Allows only a single value to be selected.
    if isempty(IncludeSSVars)
        [IncludeSSVars] = uiNameSelect({'Yes','No'},'Include SS vars?',SingleSelect);
        if strcmpi(IncludeSSVars,'Yes')
            IncludeSSVars=1;
        else
            IncludeSSVars=0;
        end
    end
elseif strcmpi(SSVarNames,'NA_EMPTY')
    IncludeSSVars=0;
    SSVarNames=[];
else
    IncludeSSVars=1;
end

if IncludeSSVars==1
    [SSVars_data_all,Use_fmriprep_table,SSVarsSubSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table,...
        'LoadVarName','SS_Vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','SSVars',...
        'AnalysisName','SSVars_All',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    Behavior_data_all=Behavior_data_all(SSVarsSubSampleInd==1,:);
    AllSSVarsNames=SSVars_data_all.Properties.VariableNames;  
    if isempty(SSVarNames)
        SSVarNames=uiNameSelect(AllSSVarsNames,'Select SS varaibles to compare.',0);
    end
    SSVars_data_all=SSVars_data_all(:,SSVarNames);
else
    SSVars_data_all=[];
end
AnalysisParameters.SSVarNames=SSVarNames;
AnalysisParameters.IncludeSSVars=IncludeSSVars;
All_Variables=[Behavior_data_all,SSVars_data_all];
All_Variables=array2table(double(table2array(All_Variables)),'VariableNames',All_Variables.Properties.VariableNames);
AnalysisParameters.All_Variables=All_Variables;
beh_var_names=All_Variables.Properties.VariableNames(:);
numvars=size(beh_var_names,1);
%% Identify parcellations availible and select ones to run analysis on.
if isempty(ParcelNames)
    AllParcelNames(ismember(AllParcelNames,'Searchlight'),:)=[];
    ParcelNames=uiNameSelect(AllParcelNames,'Select parcellations to run.');
end
AnalysisParameters.ParcelNames=ParcelNames;
numParcels=length(ParcelNames);
if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end

%% Select motion or confound variables for mixed effects model
if isempty(motion_var_names)
    if bySS==0
        motion_var_names=uiNameSelect({'MotionByRun_fd_mean','MotionByRun_fd_max','MotionByRun_rot_mean','MotionByRun_rot_max','MotionByRun_trans_mean','MotionByRun_trans_max','MotionByRun_nsso_sum'},'Select motion variables to include',0);
    else
        motion_var_names=uiNameSelect({'MotionBySS_fd_mean','MotionBySS_fd_max','MotionBySS_rot_mean','MotionBySS_rot_max','MotionBySS_trans_mean','MotionBySS_trans_max','MotionBySS_nsso_sum'},'Select motion variables to include',0);
    end
end
OnlyLMEScatter=0;
AnalysisParameters.motion_var_names=motion_var_names;
All_MotionStats=Use_fmriprep_table(:,motion_var_names);
All_MotionStats=array2table(double(table2array(All_MotionStats)),'VariableNames',All_MotionStats.Properties.VariableNames);
if isempty(RFX_Intercepts)
    RFX_Intercepts=uiNameSelect(Use_fmriprep_table.Properties.VariableNames,'Select random effect intercepts to include:',0);
elseif strcmpi(RFX_Intercepts,'NA_EMPTY')
    RFX_Intercepts=[];
end
AnalysisParameters.RFX_Intercepts=RFX_Intercepts;
if isempty(RFX_Slopes)
    RFX_Slopes=uiNameSelect(Use_fmriprep_table.Properties.VariableNames,'Select random effect slopes to include:',0);
elseif strcmpi(RFX_Slopes,'NA_EMPTY')
    RFX_Slopes=[];
end
AnalysisParameters.RFX_Slopes=RFX_Slopes;
lmem_formula=['behavVar ~ brainVar + '];
if ~isempty(motion_var_names)
    for i = 1:length(motion_var_names)
        lmem_formula =  [lmem_formula,motion_var_names{i,1},' + '];
    end
end

if ~isempty(RFX_Intercepts)
    OnlyLMEScatter=1;
    for i = 1:length(RFX_Intercepts)
        lmem_formula =  [lmem_formula,'(1|',RFX_Intercepts{i,1},')',' + '];
    end
end  
if ~isempty(RFX_Slopes)
    OnlyLMEScatter=1;
    for i = 1:length(RFX_Slopes)
        lmem_formula =  [lmem_formula,'(1 + brainVar|',RFX_Slopes{i,1},')',' + '];
    end
end  

if strcmp(lmem_formula(1,end-2:end),' + ')
    lmem_formula(:,end-2:end)=[];
end
BaseTable=[Use_fmriprep_table(:,RFX_Intercepts),Use_fmriprep_table(:,RFX_Slopes),All_MotionStats,All_Variables,Use_fmriprep_table(:,AveragingVar)];
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

    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        [NodeDataAll,~,NodeDataSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table,...
            'LoadVarName',NodeVarName,...
            'LoadVarFormat','Table',...
            'DataType','ByParcellation',...
            'AnalysisType',NodeAnalysisType,...
            'AnalysisName',NodeAnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'RowCompile',NodeRowCompileVar,...
            'ParcelName',parcelName);
        BaseTable_Nodes=BaseTable(NodeDataSampleInd==1,:);
        numNodes=length(UseLabels);
        if size(NodeDataAll,2)~=numNodes
            NodeDataAll=NodeDataAll';
            if size(NodeDataAll,2)~=numNodes
                disp('Data Alignment Error for Node compile-- please debug')
                continue
            end
        end
        if istable(NodeDataAll)
            NodeDataAll=table2array(NodeDataAll);
        end
        NodeDataAll=array2table(double(NodeDataAll),'VariableNames',UseLabels);       
        if ~isempty(RFX_Slopes)
            GroupingVar=BaseTable_Nodes.(RFX_Slopes{1,1});
        elseif ~isempty(RFX_Intercepts)
            GroupingVar=BaseTable_Nodes.(RFX_Intercepts{1,1});
        else
            GroupingVar=[];
        end
    end
    if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
        [labelPairs1]=labels2uppertriuvectorlabels(UseLabels);
        [EdgeDataAll,~,EdgeDataSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table,...
            'LoadVarName',EdgeVarName,...
            'LoadVarFormat','Table',...
            'DataType','ByParcellation',...
            'AnalysisType',EdgeAnalysisType,...
            'AnalysisName',EdgeAnalysisName,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Array',...
            'RowCompile',EdgeRowCompileVar,...            
            'ParcelName',parcelName);

        BaseTable_Edges=BaseTable(EdgeDataSampleInd==1,:);
        numEdges=length(labelPairs1);
        if ndims(EdgeDataAll)>3
            disp('Too many dimensions in Edge data-- please debug')
            continue            
        elseif ndims(EdgeDataAll)==3 && size(EdgeDataAll,1)==size(EdgeDataAll,2)
            EdgeDataAll=mat2uppertriuvectormat(EdgeDataAll);
        end     
        
        if size(EdgeDataAll,2)~=numEdges
            EdgeDataAll=EdgeDataAll';
            if size(EdgeDataAll,2)~=numEdges
                disp('Data Alignment Error for Edge compile-- please debug')
                continue
            end
        end
        if istable(EdgeDataAll)
            EdgeDataAll=table2array(EdgeDataAll);
        end        
        EdgeDataAll=array2table(double(EdgeDataAll),'VariableNames',labelPairs1);
        if ~isempty(RFX_Slopes)
            GroupingVar=BaseTable_Nodes.(RFX_Slopes{1,1});
        elseif ~isempty(RFX_Intercepts)
            GroupingVar=BaseTable_Nodes.(RFX_Intercepts{1,1});
        else
            GroupingVar=[];
        end        
    end    

    %% Run mixed effects model 4 RF
    for behavVarNum=1:numvars
        disp([parcelName,': ',beh_var_names{behavVarNum,1}])
        beh_lmem_formula=strrep(lmem_formula,'behavVar',beh_var_names{behavVarNum,1});
        if exist([GroupAnalysisDir,beh_var_names{behavVarNum,1},'.mat'],'file') && Overwrite == 0
            disp('File exists and overwrite = 0: skipping')
            continue
        end
        if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
            Node_Table=BaseTable_Nodes(:,[beh_var_names(behavVarNum,1);motion_var_names;RFX_Intercepts;RFX_Slopes;AveragingVar]);
            Temp_Results_Node=cell(length(UseLabels),1);
            Results_Node=[];
            if ~any(~isnan(table2array(BaseTable_Nodes(:,beh_var_names(behavVarNum,1)))))
                continue
            end
            if ~isempty(AveragingVar)
                
                Use_Node_Table=grpstats(Node_Table,AveragingVar{1,1});
                Use_Node_Table(:,AveragingVar)=[];
                Use_Node_Table(:,{'GroupCount'})=[];
                Use_Node_Table.Properties.VariableNames=strrepCell(Use_Node_Table.Properties.VariableNames,'mean_','');
                
                Use_NodeDataAll=[NodeDataAll,Node_Table(:,AveragingVar)];
                Use_NodeDataAll=grpstats(Use_NodeDataAll,AveragingVar{1,1});              
                Use_NodeDataAll(:,AveragingVar)=[];
                Use_NodeDataAll(:,{'GroupCount'})=[];
                Use_NodeDataAll.Properties.VariableNames=strrepCell(Use_NodeDataAll.Properties.VariableNames,'mean_','');
            else
                Use_Node_Table=Node_Table;
                Use_NodeDataAll=NodeDataAll;
            end
            % Run LMEM on RF data
            for NodeNum=1:length(UseLabels)
                Node_lmem_formula=strrep(beh_lmem_formula,'brainVar',UseLabels{NodeNum,1});
                if ~any(~isnan(table2array(Use_NodeDataAll(:,UseLabels{NodeNum,1}))))
                    continue
                end 
                try
                    [Temp_Results_Node{NodeNum,1}] = LME_PlusContinuousStats([Use_Node_Table,Use_NodeDataAll(:,UseLabels{NodeNum,1})],Node_lmem_formula);
                catch
                    Temp_Results_Node{NodeNum,1}=[];
                end
            end
            Results_Node_Confound.Intercept = [];
            Results_Node_Confound.Intercept = [];
            for i = 1:length(motion_var_names)
                Results_Node_Confound.(motion_var_names{i,1})=[];
            end
            for NodeNum=1:length(UseLabels)
                if ~isempty(Temp_Results_Node{NodeNum,1})
                    Results_Node=[Results_Node;Temp_Results_Node{NodeNum,1}(UseLabels{NodeNum,1},:)];
                    Results_Node_Confound.Intercept=[Results_Node_Confound.Intercept;Temp_Results_Node{NodeNum,1}('Intercept',:)];
                    Results_Node_Confound.Intercept.Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                    for i = 1:length(motion_var_names)
                        Results_Node_Confound.(motion_var_names{i,1})=[Results_Node_Confound.(motion_var_names{i,1});Temp_Results_Node{NodeNum,1}(motion_var_names{i,1},:)];
                        Results_Node_Confound.(motion_var_names{i,1}).Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                    end
                else
                    nanfill=array2table(nan(1,size(Results_Node,2)),'VariableNames',Results_Node.Properties.VariableNames);
                    Results_Node=[Results_Node;nanfill];
                    Results_Node.Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                    Results_Node_Confound.Intercept=[Results_Node_Confound.Intercept;nanfill];
                    Results_Node_Confound.Intercept.Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                    for i = 1:length(motion_var_names)
                        Results_Node_Confound.(motion_var_names{i,1})=[Results_Node_Confound.(motion_var_names{i,1});nanfill];
                        Results_Node_Confound.(motion_var_names{i,1}).Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                    end  
                end
            end   
            BehLabel=beh_var_names{behavVarNum,1};
            BehavVals=Use_Node_Table.(BehLabel);
        else
            Results_Node=[];
            Use_NodeDataAll=[];
            Results_Node_Confound=[];
        end
        % Run LMEM on RC data
        if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
            Edge_Table=BaseTable_Edges(:,[beh_var_names(behavVarNum,1);motion_var_names;RFX_Intercepts;RFX_Slopes;AveragingVar]);
            Temp_Results_Edge=cell(length(labelPairs1),1);
            Results_Edge=[];
             if ~isempty(AveragingVar)
                Use_Edge_Table=grpstats(Edge_Table,AveragingVar{1,1});
                Use_Edge_Table(:,AveragingVar)=[];
                Use_Edge_Table(:,{'GroupCount'})=[];
                Use_Edge_Table.Properties.VariableNames=strrepCell(Use_Edge_Table.Properties.VariableNames,'mean_','');
                
                Use_EdgeDataAll=[EdgeDataAll,Edge_Table(:,AveragingVar)];
                Use_EdgeDataAll=grpstats(Use_EdgeDataAll,AveragingVar{1,1});              
                Use_EdgeDataAll(:,AveragingVar)=[];
                Use_EdgeDataAll(:,{'GroupCount'})=[];
                Use_EdgeDataAll.Properties.VariableNames=strrepCell(Use_EdgeDataAll.Properties.VariableNames,'mean_','');
            else
                Use_Edge_Table=Edge_Table;
                Use_EdgeDataAll=EdgeDataAll;
            end           
            parfor EdgeNum=1:length(labelPairs1)
                Edge_lmem_formula=strrep(beh_lmem_formula,'brainVar',labelPairs1{EdgeNum,1});
                if ~any(~isnan(table2array(Use_EdgeDataAll(:,labelPairs1{EdgeNum,1}))))
                    continue
                end    
                try
                    [Temp_Results_Edge{EdgeNum,1}] = LME_PlusContinuousStats([Use_Edge_Table,Use_EdgeDataAll(:,labelPairs1{EdgeNum,1})],Edge_lmem_formula);  
                catch
                    Temp_Results_Edge{EdgeNum,1}=[];
                end    
            end 
            Results_Edge_Confound.Intercept = [];
            for i = 1:length(motion_var_names)
                Results_Edge_Confound.(motion_var_names{i,1})=[];
            end  
             for EdgeNum=1:length(labelPairs1)
                if ~isempty(Temp_Results_Edge{EdgeNum,1})                
                    Results_Edge=[Results_Edge;Temp_Results_Edge{EdgeNum,1}(labelPairs1{EdgeNum,1},:)];
                    Results_Edge_Confound.Intercept=[Results_Edge_Confound.Intercept;Temp_Results_Edge{EdgeNum,1}('Intercept',:)];
                    Results_Edge_Confound.Intercept.Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                    for i = 1:length(motion_var_names)
                        Results_Edge_Confound.(motion_var_names{i,1})=[Results_Edge_Confound.(motion_var_names{i,1});Temp_Results_Edge{EdgeNum,1}(motion_var_names{i,1},:)];
                        Results_Edge_Confound.(motion_var_names{i,1}).Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                    end
                else
                    nanfill=array2table(nan(1,size(Results_Edge,2)),'VariableNames',Results_Edge.Properties.VariableNames);
                    Results_Edge=[Results_Edge;nanfill];
                    Results_Edge.Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                    Results_Edge_Confound.Intercept=[Results_Edge_Confound.Intercept;nanfill];
                    Results_Edge_Confound.Intercept.Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                    for i = 1:length(motion_var_names)
                        Results_Edge_Confound.(motion_var_names{i,1})=[Results_Edge_Confound.(motion_var_names{i,1});nanfill];
                        Results_Edge_Confound.(motion_var_names{i,1}).Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                    end  
                end
             end 
             BehLabel=beh_var_names{behavVarNum,1};
             BehavVals=Use_Edge_Table.(BehLabel);
         
        else
            Results_Edge=[];
            Use_EdgeDataAll=[];
            Results_Edge_Confound=[];
        end
        if MakeFigs == 1
            [Edge_Degree_Table] = IndvDiffParcellationSummaryFigures(Results_Node,Results_Edge,Use_NodeDataAll,Use_EdgeDataAll,BehavVals,BehLabel,NodeVarName,EdgeVarName,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'GroupingVar',GroupingVar,'OnlyLMEScatter',OnlyLMEScatter);
        else
            Edge_Degree_Table=[];
        end
        if ~exist(GroupAnalysisDir,'file')
            mkdir(GroupAnalysisDir);
        end
        save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
        save([GroupAnalysisDir,beh_var_names{behavVarNum,1}],'Results_Node','Results_Node_Confound','Results_Edge','Results_Edge_Confound','Edge_Degree_Table');
        toc       
    end  
end
