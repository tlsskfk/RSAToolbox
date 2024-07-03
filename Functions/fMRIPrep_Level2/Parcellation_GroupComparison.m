function [AnalysisParameters] = Parcellation_GroupComparison(ExperimentsDir,varargin)
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
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    SSVarNames=AnalysisParameters.SSVarNames;
    IncludeSSVars=AnalysisParameters.IncludeSSVars;    
    NodeAnalysisType=AnalysisParameters.NodeAnalysisType;
    NodeAnalysisName=AnalysisParameters.NodeAnalysisName;
    EdgeAnalysisType=AnalysisParameters.EdgeAnalysisType;
    EdgeAnalysisName=AnalysisParameters.EdgeAnalysisName;
    beh_var_names=AnalysisParameters.beh_var_names;
    motion_var_names=AnalysisParameters.motion_var_names;
    fmriprep_table_names=AnalysisParameters.fmriprep_table_names;
%    fmriprep_tables=AnalysisParameters.fmriprep_tables;
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

if isempty(fmriprep_table_names)
    NumGroups=uiEnterName('2','Enter # of groups to compare');
    NumGroups=str2num(NumGroups);
    for i = 1:NumGroups
        [~,~,fmriprep_table_names{i,1}] = load_fmriprep_table;
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
% fmriprep_tables=AnalysisParameters.fmriprep_tables;
GroupVar=array2table(GroupVar,'VariableNames',{'GroupVar'});
fmriprep_table=[fmriprep_tables,GroupVar];

AnalysisParameters.fmriprep_table_names=fmriprep_table_names;
AnalysisParameters.fmriprep_table=fmriprep_table;
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupDiff_GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupDiff_GroupAnalysis/','/GroupDiff_GroupFigures/');
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

if isempty(RunVarTypes)
    SingleSelect=1; %Allows only a single value to be selected.
    [RunVarTypes] = uiNameSelect({'Nodes','Edges','Both'},'GroupDiff measures on parcel...',SingleSelect);
end
AnalysisParameters.RunVarTypes=RunVarTypes;
if isempty(NodeVarName)
   NodeVarName=cell(NumGroups,1);
end
if isempty(EdgeVarName)
   EdgeVarName=cell(NumGroups,1);
end
if isempty(NodeAnalysisType)
   NodeAnalysisType=cell(NumGroups,1);
end
if isempty(EdgeAnalysisType)
   EdgeAnalysisType=cell(NumGroups,1);
end
if isempty(NodeAnalysisName)
   NodeAnalysisName=cell(NumGroups,1);
end
if isempty(EdgeAnalysisName)
   EdgeAnalysisName=cell(NumGroups,1);
end
if isempty(NodeRowCompileVar)
   NodeRowCompileVar=cell(NumGroups,1);
end
if isempty(EdgeRowCompileVar)
   EdgeRowCompileVar=cell(NumGroups,1);
end
filePaths_NodeData=cell(NumGroups,1);
filePaths_EdgeData=cell(NumGroups,1);
%% Compile filepaths for input files for the analysis
for i = 1:NumGroups
    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        [filePaths_NodeData{i,1},NodeAnalysisType{i,1},NodeAnalysisName{i,1}] = BIDsDirSearch(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),'SubjectOrRun',SubjectOrRun,'AnalysisType',NodeAnalysisType{i,1},'AnalysisName',NodeAnalysisName{i,1},'TitleTextName',['Select source of node values for',newline,fmriprep_table_names{i,2}]);
        AllParcelNamesNodes=filePaths_NodeData{i,1}.Properties.VariableNames;
        if isempty(NodeVarName{i,1})
            tempNodeVarName=[];
            for j = 1:size(filePaths_NodeData{i,1},1)
                try
                    variableInfo = who('-file', filePaths_NodeData{i,1}.(AllParcelNamesNodes{1,1}){j,1});
                    tempNodeVarName=unique([tempNodeVarName;variableInfo]);
                catch
                    continue
                end
            end
            [NodeVarName{i,1}] = uiNameSelect(tempNodeVarName,['Select node variable for',newline,fmriprep_table_names{i,2}] ,1);    
        end
        if isempty(NodeRowCompileVar{i,1})
            for j = 1:size(filePaths_NodeData{i,1},1)
                try
                    tempNodeData=load(filePaths_NodeData{i,1}.(AllParcelNamesNodes{1,1}){j,1},NodeVarName{i,1});
                    tempNodeData=tempNodeData.(NodeVarName{i,1});
                    if istable(tempNodeData)
                        if ~isempty(tempNodeData.Properties.RowNames)
                            NodeRowCompileVar{i,1}=uiNameSelect(tempNodeData.Properties.RowNames,['Select node row variable for',newline,fmriprep_table_names{i,2}],1);
                            break
                        end
                    end    
                catch
                    continue
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
    if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
        [filePaths_EdgeData{i,1},EdgeAnalysisType{i,1},EdgeAnalysisName{i,1}] = BIDsDirSearch(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),'SubjectOrRun',SubjectOrRun,'AnalysisType',EdgeAnalysisType{i,1},'AnalysisName',EdgeAnalysisName{i,1},'TitleTextName',['Select source of edge values for',newline,fmriprep_table_names{i,2}]);
        AllParcelNamesEdges=filePaths_EdgeData{i,1}.Properties.VariableNames;
        if isempty(EdgeVarName{i,1})
            tempEdgeVarName=[];
            for j = 1:size(filePaths_EdgeData{i,1},1)
                try
                    variableInfo = who('-file', filePaths_EdgeData{i,1}.(AllParcelNamesEdges{1,1}){j,1});
                    tempEdgeVarName=unique([tempEdgeVarName;variableInfo]);
                catch
                    continue
                end
            end
            [EdgeVarName{i,1}] = uiNameSelect(tempEdgeVarName,['Select edge variable for',newline,fmriprep_table_names{i,2}] ,1);
        end
        if isempty(EdgeRowCompileVar{i,1})
            for j = 1:size(filePaths_EdgeData{i,1},1)
                try
                    tempEdgeData=load(filePaths_EdgeData{i,1}.(AllParcelNamesEdges{1,1}){j,1},EdgeVarName{i,1});
                    tempEdgeData=tempEdgeData.(EdgeVarName{i,1});
                    if istable(tempEdgeData)                    
                        if ~isempty(tempEdgeData.Properties.RowNames)
                            if ~any(~ismember(tempEdgeData.Properties.RowNames,tempEdgeData.Properties.VariableNames)) && length(tempEdgeData.Properties.RowNames)==length(tempEdgeData.Properties.VariableNames)
                                break
                            end
                            EdgeRowCompileVar{i,1}=uiNameSelect(tempEdgeData.Properties.RowNames,['Select Edge row variable for',newline,fmriprep_table_names{i,2}],1);
                            break
                        end
                    end    
                catch
                    continue
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
end
AnalysisParameters.NodeAnalysisType=NodeAnalysisType;
AnalysisParameters.NodeAnalysisName=NodeAnalysisName;
AnalysisParameters.NodeVarName=NodeVarName;   
AnalysisParameters.NodeRowCompileVar=NodeRowCompileVar;   
AnalysisParameters.EdgeAnalysisType=EdgeAnalysisType;
AnalysisParameters.EdgeAnalysisName=EdgeAnalysisName;
AnalysisParameters.EdgeVarName=EdgeVarName; 
AnalysisParameters.EdgeRowCompileVar=EdgeRowCompileVar; 

if strcmpi(RunVarTypes,'Nodes')
    AllParcelNames=AllParcelNamesNodes;
    AnalysisType=['GroupDiff_',NodeAnalysisType{1,1}];
    AnalysisNamesSuffix=[NodeAnalysisName{1,1},'_',NodeVarName{1,1},NodeRowCompileVar{1,1}];
elseif strcmpi(RunVarTypes,'Edges')
    AllParcelNames=AllParcelNamesEdges;
    AnalysisType=['GroupDiff_',EdgeAnalysisType{1,1}];
    AnalysisNamesSuffix=[EdgeAnalysisName{1,1},'_',EdgeVarName{1,1},EdgeRowCompileVar{1,1}];   
elseif strcmpi(RunVarTypes,'Both')
    AllParcelNames=AllParcelNamesNodes(ismember(AllParcelNamesNodes,AllParcelNamesEdges)==1);
    AnalysisType=['GroupDiff_',NodeAnalysisType{1,1},'_',EdgeAnalysisType{1,1}];
    if strcmpi(NodeAnalysisName{1,1},EdgeAnalysisName{1,1})
        AnalysisNamesSuffix=[NodeAnalysisName{1,1},'_',NodeVarName{1,1},NodeRowCompileVar{1,1},'_',EdgeVarName{1,1},EdgeRowCompileVar{1,1}];
    else
        AnalysisNamesSuffix=[NodeAnalysisName{1,1},'_',EdgeAnalysisName{1,1},'_',NodeVarName{1,1},NodeRowCompileVar{1,1},'_',EdgeVarName{1,1},EdgeRowCompileVar{1,1}];
    end
end
if DefaultName==0
    if isempty(AnalysisName)
        AnalysisName=uiEnterName(['GroupDiff_',AnalysisNamesSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    end
    ComparisonNames=uiEnterName([ComparisonNames],['Enter name for ',AnalysisType,newline,'analysis below:']);
else
    AnalysisName=['GroupDiff_',AnalysisNamesSuffix];
end

%%Compile Behavioral Data
Use_fmriprep_table=[];
Behavior_data_all=[];

for i = 1:NumGroups
    try
        try
        [tempBehavior_data_all,tempUse_fmriprep_table,~] = CompileND(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),...
            'LoadVarName','beh_vars',...
            'LoadVarFormat','Table',...
            'DataType','Other',...
            'AnalysisType','beh_vars',...
            'AnalysisName','Overall',...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Table',...
            'ParcelName','None');
        catch
        [tempBehavior_data_all,tempUse_fmriprep_table,~] = CompileND(ExperimentsDir,fmriprep_table(fmriprep_table.GroupVar==i,:),...
            'LoadVarName','SS_Vars',...
            'LoadVarFormat','Table',...
            'DataType','Other',...
            'AnalysisType','SSVars',...
            'AnalysisName','SSVars_All',...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Table',...
            'ParcelName','None');
        end
        Use_fmriprep_table=[Use_fmriprep_table;tempUse_fmriprep_table];
        Behavior_data_all=[Behavior_data_all;tempBehavior_data_all];    
    catch
        Use_fmriprep_table=[Use_fmriprep_table;fmriprep_table.GroupVar==i];
    end
end

if ~isempty(Behavior_data_all)
All_beh_var_names=Behavior_data_all.Properties.VariableNames(:);
if isempty(beh_var_names)
    beh_var_names=uiNameSelect([{'none'};All_beh_var_names],'Select behavioral varaibles for covariates.',0);
end
end
if ~isempty(beh_var_names)
    if ~iscell(beh_var_names)
        beh_var_names={beh_var_names};
    end
    if any(strcmpi(beh_var_names,'none')) || any(strcmpi(beh_var_names,'NA_EMPTY'))
        beh_var_names=[];
    end
else
    beh_var_names=[];
end
    
    
AnalysisParameters.beh_var_names=beh_var_names;
Behavior_data_all=Behavior_data_all(:,beh_var_names);

if isempty(SSVarNames)
    SingleSelect=1; %Allows only a single value to be selected.
    if isempty(IncludeSSVars)
        [IncludeSSVars] = uiNameSelect({'Yes','No'},'Include SS vars as covariates?',SingleSelect);
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
    Use_fmriprep_table2=Use_fmriprep_table;
    Use_fmriprep_table=[];
    SSVarsSubSampleInd=[];
    SSVars_data_all=[];
    for i = 1:NumGroups
        [tempSSVars_data_all,tempUse_fmriprep_table,tempSSVarsSubSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table2(Use_fmriprep_table2.GroupVar==i,:),...
            'LoadVarName','SS_Vars',...
            'LoadVarFormat','Table',...
            'DataType','Other',...
            'AnalysisType','SSVars',...
            'AnalysisName','SSVars_All',...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray','Table',...
            'ParcelName','None');
        Use_fmriprep_table=[Use_fmriprep_table;tempUse_fmriprep_table];
        SSVarsSubSampleInd=[SSVarsSubSampleInd;tempSSVarsSubSampleInd];
        SSVars_data_all=[SSVars_data_all;tempSSVars_data_all];
    end
    Behavior_data_all=Behavior_data_all(SSVarsSubSampleInd==1,:);
    AllSSVarsNames=SSVars_data_all.Properties.VariableNames;  
    if isempty(SSVarNames)
        SSVarNames=uiNameSelect(AllSSVarsNames,'Select SS varaibles for covariates.',0);
    end
    SSVars_data_all=SSVars_data_all(:,SSVarNames);
else
    SSVars_data_all=[];
end
AnalysisParameters.SSVarNames=SSVarNames;
AnalysisParameters.IncludeSSVars=IncludeSSVars;
try
    All_Variables=[Behavior_data_all,SSVars_data_all];
    All_Variables=array2table(double(table2array(All_Variables)),'VariableNames',All_Variables.Properties.VariableNames);
    AnalysisParameters.All_Variables=All_Variables;
    beh_var_names=All_Variables.Properties.VariableNames(:);
    numvars=size(beh_var_names,1);
catch
    All_Variables=[];
end
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
        motion_var_names=uiNameSelect({'none','MotionByRun_fd_mean','MotionByRun_fd_max','MotionByRun_rot_mean','MotionByRun_rot_max','MotionByRun_trans_mean','MotionByRun_trans_max','MotionByRun_nsso_sum'},'Select motion variables to include',0);
    else
        motion_var_names=uiNameSelect({'none','MotionBySS_fd_mean','MotionBySS_fd_max','MotionBySS_rot_mean','MotionBySS_rot_max','MotionBySS_trans_mean','MotionBySS_trans_max','MotionBySS_nsso_sum'},'Select motion variables to include',0);
    end
end
if any(strcmpi(motion_var_names,'none')) || any(strcmpi(motion_var_names,'NA_EMPTY'))
    motion_var_names=[];
end
OnlyLMEScatter=0;
AnalysisParameters.motion_var_names=motion_var_names;
if ~isempty(motion_var_names)
All_MotionStats=Use_fmriprep_table(:,motion_var_names);
All_MotionStats=array2table(double(table2array(All_MotionStats)),'VariableNames',All_MotionStats.Properties.VariableNames);
else
    All_MotionStats=[];
end
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
lmem_formula=['brainVar ~ GroupVar + '];
if ~isempty(beh_var_names)
    for i = 1:length(beh_var_names)
        lmem_formula =  [lmem_formula,beh_var_names{i,1},' + '];
    end
end
if ~isempty(SSVarNames)
    for i = 1:length(SSVarNames)
        lmem_formula =  [lmem_formula,SSVarNames{i,1},' + '];
    end
end
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
        lmem_formula =  [lmem_formula,'(1 + GroupVar|',RFX_Slopes{i,1},')',' + '];
    end
end  

if strcmp(lmem_formula(1,end-2:end),' + ')
    lmem_formula(:,end-2:end)=[];
end

BaseTable=[Use_fmriprep_table(:,{'GroupVar'}),Use_fmriprep_table(:,RFX_Intercepts),Use_fmriprep_table(:,RFX_Slopes)];
if ~isempty(All_MotionStats)
    BaseTable=[BaseTable,All_MotionStats];
end
if ~isempty(All_MotionStats)
    BaseTable=[BaseTable,All_Variables];
end
for parcelNum=1:numParcels
    tic
    parcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,ComparisonNames,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
        GroupFiguresDir=[FigureDir,ComparisonNames,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,ComparisonNames,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,ComparisonNames,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
        GroupFiguresDir=[FigureDir,ComparisonNames,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,ComparisonNames,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
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
        filePaths=BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisName',parcelName,'AnalysisType','CoordParcels','SubjectOrRun','Subject');
        filePaths=filePaths.(parcelName);
        for loadCt=1:length(filePaths)
            try
                load(filePaths{loadCt,1},'UseMask','UseLabels')
            catch
                continue
            end
            break
        end
    end
    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        NodeDataAll=[];
        NodeDataSampleInd=[];
        for i = 1:NumGroups
            [tempNodeDataAll,~,tempNodeDataSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table(Use_fmriprep_table.GroupVar==i,:),...
                'LoadVarName',NodeVarName{i,1},...
                'LoadVarFormat','Table',...
                'DataType','ByParcellation',...
                'AnalysisType',NodeAnalysisType{i,1},...
                'AnalysisName',NodeAnalysisName{i,1},...
                'SubjectOrRun',SubjectOrRun,...
                'TableOrArray','Array',...
                'RowCompile',NodeRowCompileVar{i,1},...
                'ParcelName',parcelName);           
            NodeDataAll=[NodeDataAll;tempNodeDataAll];
            NodeDataSampleInd=[NodeDataSampleInd;tempNodeDataSampleInd];            
        end    
          
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
        numEdges=length(labelPairs1);
        EdgeDataAll=[];
        EdgeDataSampleInd=[];        
        for i = 1:NumGroups
            [tempEdgeDataAll,~,tempEdgeDataSampleInd] = CompileND(ExperimentsDir,Use_fmriprep_table(Use_fmriprep_table.GroupVar==i,:),...
                'LoadVarName',EdgeVarName{i,1},...
                'LoadVarFormat','Table',...
                'DataType','ByParcellation',...
                'AnalysisType',EdgeAnalysisType{i,1},...
                'AnalysisName',EdgeAnalysisName{i,1},...
                'SubjectOrRun',SubjectOrRun,...
                'TableOrArray','Array',...
                'RowCompile',EdgeRowCompileVar{i,1},...            
                'ParcelName',parcelName);
            if ndims(tempEdgeDataAll)==3
                tempEdgeDataAll=mat2uppertriuvectormat(tempEdgeDataAll);
            elseif ndims(tempEdgeDataAll)>3
                disp('Too many dimensions in Edge data-- please debug')
                continue
            end
            if size(tempEdgeDataAll,2)~=numEdges
                tempEdgeDataAll=tempEdgeDataAll';
                if size(tempEdgeDataAll,2)~=numEdges
                    disp('Data Alignment Error for Edge compile-- please debug')
                    continue
                end
            end            
            EdgeDataAll=[EdgeDataAll;tempEdgeDataAll];
            EdgeDataSampleInd=[EdgeDataSampleInd;tempEdgeDataSampleInd];            
        end
        BaseTable_Edges=BaseTable(EdgeDataSampleInd==1,:);
        
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
    GroupLabels=fmriprep_table_names(:,2);
    [~,groupPairLabels]=labels2uppertriuvectorlabels(GroupLabels);
    groupPairLabels=strrepCell(groupPairLabels,'_2_','vs');
    %% Run mixed effects model 4 RF
    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        % Run LMEM on RF data
        Temp_Results_Node=cell(length(UseLabels),1);
        Temp_PairwisePermTable=cell(length(UseLabels),1);
        Temp_IndvPermTable=cell(length(UseLabels),1);
        Results_Node=[];
        Perm_Pairwise_Results_Node=struct;
        parfor NodeNum=1:length(UseLabels)
            Node_lmem_formula=strrep(lmem_formula,'brainVar',UseLabels{NodeNum,1});
            if ~any(~isnan(table2array(NodeDataAll(:,UseLabels{NodeNum,1}))))
                continue
            end 
            try
                [Temp_Results_Node{NodeNum,1},Temp_PairwisePermTable{NodeNum,1},Temp_IndvPermTable{NodeNum,1}] = LME_PlusGroupStats([BaseTable_Nodes,NodeDataAll(:,UseLabels{NodeNum,1})],Node_lmem_formula,GroupLabels);
            catch
                Temp_Results_Node{NodeNum,1}=[];
                Temp_PairwisePermTable{NodeNum,1}=[];
                Temp_IndvPermTable{NodeNum,1}=[];
            end
        end
        Results_Node_Confound.Intercept = [];
        Results_Node_Confound.Intercept = [];
        for i = 1:length(motion_var_names)
            Results_Node_Confound.(motion_var_names{i,1})=[];
        end
        for i = 1:length(groupPairLabels)
            Perm_Pairwise_Results_Node.(groupPairLabels{i,1})=[];
        end   
        for i = 1:length(GroupLabels)
            Perm_Indv_Results_Node.(GroupLabels{i,1})=[];
        end         
        for NodeNum=1:length(UseLabels)
            if ~isempty(Temp_Results_Node{NodeNum,1})
                Temp_Results_Node{NodeNum,1}.Properties.RowNames=strrepCell(Temp_Results_Node{NodeNum,1}.Properties.RowNames,'GroupVar',UseLabels{NodeNum,1});
                Results_Node=[Results_Node;Temp_Results_Node{NodeNum,1}(UseLabels{NodeNum,1},:)];
                Results_Node_Confound.Intercept=[Results_Node_Confound.Intercept;Temp_Results_Node{NodeNum,1}('Intercept',:)];
                Results_Node_Confound.Intercept.Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                for i = 1:length(motion_var_names)
                    Results_Node_Confound.(motion_var_names{i,1})=[Results_Node_Confound.(motion_var_names{i,1});Temp_Results_Node{NodeNum,1}(motion_var_names{i,1},:)];
                    Results_Node_Confound.(motion_var_names{i,1}).Properties.RowNames{NodeNum,1}=UseLabels{NodeNum,1};
                end
                for i = 1:length(groupPairLabels)
                    Perm_Pairwise_Results_Node.(groupPairLabels{i,1})=[Perm_Pairwise_Results_Node.(groupPairLabels{i,1});Temp_PairwisePermTable{NodeNum,1}(groupPairLabels{i,1},:)];
                    Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).Properties.RowNames=strrepCell(Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).Properties.RowNames,groupPairLabels{i,1},UseLabels{NodeNum,1});
                end
                for i = 1:length(GroupLabels)
                    Perm_Indv_Results_Node.(GroupLabels{i,1})=[Perm_Indv_Results_Node.(GroupLabels{i,1});Temp_IndvPermTable{NodeNum,1}(GroupLabels{i,1},:)];  
                    Perm_Indv_Results_Node.(GroupLabels{i,1}).Properties.RowNames=strrepCell(Perm_Indv_Results_Node.(GroupLabels{i,1}).Properties.RowNames,GroupLabels{i,1},UseLabels{NodeNum,1});                    
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
    else
        Results_Node=[];
        NodeDataAll=[];
        Results_Node_Confound=[];
    end
    % Run LMEM on RC data
    if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
        Temp_Results_Edge=cell(length(labelPairs1),1);
        Temp_PairwisePermTable=cell(length(labelPairs1),1);
        Temp_IndvPermTable=cell(length(labelPairs1),1);
        Perm_Pairwise_Results_Edge=struct;        
        Results_Edge=[];
        parfor EdgeNum=1:length(labelPairs1)
            Edge_lmem_formula=strrep(lmem_formula,'brainVar',labelPairs1{EdgeNum,1});
            if ~any(~isnan(table2array(EdgeDataAll(:,labelPairs1{EdgeNum,1}))))
                continue
            end    
            try
                [Temp_Results_Edge{EdgeNum,1},Temp_PairwisePermTable{EdgeNum,1},Temp_IndvPermTable{EdgeNum,1}] = LME_PlusGroupStats([BaseTable_Edges,EdgeDataAll(:,labelPairs1{EdgeNum,1})],Edge_lmem_formula,GroupLabels);
            catch
                Temp_Results_Edge{EdgeNum,1}=[];
                Temp_PairwisePermTable{EdgeNum,1}=[];
                Temp_IndvPermTable{EdgeNum,1}=[];               
            end    
        end 
        Results_Edge_Confound.Intercept = [];
        for i = 1:length(motion_var_names)
            Results_Edge_Confound.(motion_var_names{i,1})=[];
        end  
        for i = 1:length(groupPairLabels)
            Perm_Pairwise_Results_Edge.(groupPairLabels{i,1})=[];
        end   
        for i = 1:length(GroupLabels)
            Perm_Indv_Results_Edge.(GroupLabels{i,1})=[];
        end     
        for EdgeNum=1:length(labelPairs1)
            if ~isempty(Temp_Results_Edge{EdgeNum,1})
                Temp_Results_Edge{EdgeNum,1}.Properties.RowNames=strrepCell(Temp_Results_Edge{EdgeNum,1}.Properties.RowNames,'GroupVar',labelPairs1{EdgeNum,1});
                Results_Edge=[Results_Edge;Temp_Results_Edge{EdgeNum,1}(labelPairs1{EdgeNum,1},:)];
                Results_Edge_Confound.Intercept=[Results_Edge_Confound.Intercept;Temp_Results_Edge{EdgeNum,1}('Intercept',:)];
                Results_Edge_Confound.Intercept.Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                for i = 1:length(motion_var_names)
                    Results_Edge_Confound.(motion_var_names{i,1})=[Results_Edge_Confound.(motion_var_names{i,1});Temp_Results_Edge{EdgeNum,1}(motion_var_names{i,1},:)];
                    Results_Edge_Confound.(motion_var_names{i,1}).Properties.RowNames{EdgeNum,1}=labelPairs1{EdgeNum,1};
                end
                for i = 1:length(groupPairLabels)
                    Perm_Pairwise_Results_Edge.(groupPairLabels{i,1})=[Perm_Pairwise_Results_Edge.(groupPairLabels{i,1});Temp_PairwisePermTable{EdgeNum,1}(groupPairLabels{i,1},:)];
                    Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).Properties.RowNames=strrepCell(Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).Properties.RowNames,groupPairLabels{i,1},labelPairs1{EdgeNum,1});
                end
                for i = 1:length(GroupLabels)
                    Perm_Indv_Results_Edge.(GroupLabels{i,1})=[Perm_Indv_Results_Edge.(GroupLabels{i,1});Temp_IndvPermTable{EdgeNum,1}(GroupLabels{i,1},:)];  
                    Perm_Indv_Results_Edge.(GroupLabels{i,1}).Properties.RowNames=strrepCell(Perm_Indv_Results_Edge.(GroupLabels{i,1}).Properties.RowNames,GroupLabels{i,1},labelPairs1{EdgeNum,1});                    
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
    else
        Results_Edge=[];
        EdgeDataAll=[];
        Results_Edge_Confound=[];
    end
    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        inData_Node=cell(length(GroupLabels),length(UseLabels));
        NodeStats.Indv_tP=zeros(length(GroupLabels),length(UseLabels));
        NodeStats.Indv_permP=zeros(length(GroupLabels),length(UseLabels));
        NodeStats.Indv_val=zeros(length(GroupLabels),length(UseLabels)); 
        NodeStats.Indv_Tval=zeros(length(GroupLabels),length(UseLabels));
        for i = 1:length(GroupLabels)
            if length(UseLabels)==1
                inData_Node{i,1}=table2array(NodeDataAll(BaseTable_Nodes.GroupVar==i,:));
            else
                inData_Node(i,:)=num2cell(table2array(NodeDataAll(BaseTable_Nodes.GroupVar==i,:)),[1,length(UseLabels)]);
            end
            NodeStats.Indv_tP(i,:)=Perm_Indv_Results_Node.(GroupLabels{i,1}).Tp';
            NodeStats.Indv_permP(i,:)=Perm_Indv_Results_Node.(GroupLabels{i,1}).PermP';
            NodeStats.Indv_val(i,:)=Perm_Indv_Results_Node.(GroupLabels{i,1}).Diff';  
            NodeStats.Indv_Tval(i,:)=Perm_Indv_Results_Node.(GroupLabels{i,1}).Tval';
        end
        NodeStats.Pairwise_lmeP=zeros(length(groupPairLabels),length(UseLabels));
        NodeStats.Pairwise_tP=zeros(length(groupPairLabels),length(UseLabels));
        NodeStats.Pairwise_permP=zeros(length(groupPairLabels),length(UseLabels));
        NodeStats.Pairwise_val=zeros(length(groupPairLabels),length(UseLabels)); 
        NodeStats.Pairwise_Tval=zeros(length(groupPairLabels),length(UseLabels));  
        for i = 1:length(groupPairLabels)           
            NodeStats.Pairwise_lmeP(i,:)=Results_Node.lme_pValue';
            NodeStats.Pairwise_tP(i,:)=Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).Tp';
            NodeStats.Pairwise_permP(i,:)=Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).PermP';
            NodeStats.Pairwise_val(i,:)=Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).Diff';
            NodeStats.Pairwise_Tval(i,:)=Perm_Pairwise_Results_Node.(groupPairLabels{i,1}).Tval';
        end
    else
        inData_Node=[];
        NodeStats=[];
        NodeVarName=[{''}];
        Perm_Indv_Results_Node=[];
        Perm_Pairwise_Results_Node=[];
        Results_Node=[];
    end    
        
    if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
        inData_Edge=cell(length(GroupLabels),length(labelPairs1));
        EdgeStats.Indv_tP=zeros(length(GroupLabels),length(labelPairs1));
        EdgeStats.Indv_permP=zeros(length(GroupLabels),length(labelPairs1));
        EdgeStats.Indv_val=zeros(length(GroupLabels),length(labelPairs1));  
        EdgeStats.Indv_Tval=zeros(length(GroupLabels),length(labelPairs1));   
        for i = 1:length(GroupLabels)
            if length(UseLabels)==1
                continue
            else    
                inData_Edge(i,:)=num2cell(table2array(EdgeDataAll(BaseTable_Edges.GroupVar==i,:)),[1,length(labelPairs1)]);
            end
            EdgeStats.Indv_tP(i,:)=Perm_Indv_Results_Edge.(GroupLabels{i,1}).Tp';
            EdgeStats.Indv_permP(i,:)=Perm_Indv_Results_Edge.(GroupLabels{i,1}).PermP';
            EdgeStats.Indv_val(i,:)=Perm_Indv_Results_Edge.(GroupLabels{i,1}).Diff'; 
            EdgeStats.Indv_Tval(i,:)=Perm_Indv_Results_Edge.(GroupLabels{i,1}).Tval'; 
        end
        EdgeStats.Pairwise_lmeP=zeros(length(groupPairLabels),length(labelPairs1));
        EdgeStats.Pairwise_tP=zeros(length(groupPairLabels),length(labelPairs1));
        EdgeStats.Pairwise_permP=zeros(length(groupPairLabels),length(labelPairs1));
        EdgeStats.Pairwise_val=zeros(length(groupPairLabels),length(labelPairs1));  
        EdgeStats.Pairwise_Tval=zeros(length(groupPairLabels),length(labelPairs1)); 
        for i = 1:length(groupPairLabels)           
            EdgeStats.Pairwise_lmeP(i,:)=Results_Edge.lme_pValue';
            EdgeStats.Pairwise_tP(i,:)=Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).Tp';
            EdgeStats.Pairwise_permP(i,:)=Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).PermP';
            EdgeStats.Pairwise_val(i,:)=Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).Diff';
            EdgeStats.Pairwise_Tval(i,:)=Perm_Pairwise_Results_Edge.(groupPairLabels{i,1}).Tval';
        end
    else
        inData_Edge=[];
        EdgeStats=[];
        EdgeVarName=[{''}];
        Perm_Indv_Results_Edge=[];
        Perm_Pairwise_Results_Edge=[];        
        EdgeVarName=[{''}];
    end
    
    
    [Edge_Degree_Table,ReliabilityResults] = GroupDiffParcellationSummaryFigures(inData_Node,inData_Edge,NodeStats,EdgeStats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,GroupLabels,groupPairLabels,'vecName',NodeVarName{1,1},'matName',EdgeVarName{1,1});

    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);
    end
    save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
    save([GroupAnalysisDir,ComparisonNames],'Results_Node','Results_Node_Confound','Results_Edge','Results_Edge_Confound','Perm_Pairwise_Results_Node','Perm_Indv_Results_Node','Perm_Pairwise_Results_Edge','Perm_Indv_Results_Edge','Edge_Degree_Table','ReliabilityResults');
    toc       
end  
end
