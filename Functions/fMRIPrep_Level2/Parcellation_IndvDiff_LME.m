function [AnalysisParameters] = Parcellation_IndvDiff_LME(varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Cell containing list of table names
[ExperimentsDir] = VariableSetter('ExperimentsDir',[],varargin);
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
[featureNames] = VariableSetter('featureNames',[],varargin);
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
[RFX_Slopes] = VariableSetter('RFX_Slopes',[],varargin);
[RFX_Intercepts] = VariableSetter('RFX_Intercepts',[],varargin);
[IncludeTaskVars] = VariableSetter('IncludeTaskVars',[],varargin); %1 = yes 
[TaskVarNames] = VariableSetter('TaskVarNames',[],varargin);

[DefaultName] = VariableSetter('DefaultName',0,varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);


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

%% Input number of groups and identify corresponding directories and tablestables
if isempty(fmriprep_table_names)
    NumGroups=uiEnterName('2','Enter # of groups to include');
    NumGroups=str2num(NumGroups);
    ExperimentsDirs=cell(NumGroups,1);
    fmriprep_table_names=cell(NumGroups,2);
    subIDNames=cell(NumGroups,1);
    for i = 1:NumGroups
        [ExperimentsDirs{i,1},tempTable,fmriprep_table_names{i,1}] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
        fmriprep_table_names{i,2}=uiEnterName(fmriprep_table_names{i,1},'Enter groupName');
        SingleSelect=1; %Allows only a single value to be selected.
        [subIDNames{i,1}] = uiNameSelect(tempTable.Properties.VariableNames,'Select variable for SS names:',SingleSelect);       
    end
end

%% Load ss tables for each group and define GroupingVar
fmriprep_tables=cell(NumGroups,1);
GroupingVar=[];
ComparisonNames=[];
for i = 1:NumGroups
    [~,fmriprep_tables{i,1}] = load_fmriprep_table('ExperimentsDir',ExperimentsDirs{i,1},'fmriprep_table_name',fmriprep_table_names{i,1});
    GroupingVar=[GroupingVar;repmat(fmriprep_table_names(i,2),[height(fmriprep_tables{i,1}),1])];
    ComparisonNames=[ComparisonNames,fmriprep_table_names{i,2}];
end


%% Set features to load
if isempty(featureNames)
    NumFeatures=uiEnterName('2','Enter # of features to include');
    NumFeatures=str2num(NumFeatures);
    featureNames=cell(1,NumFeatures);
    for i = 1:NumFeatures
        featureNames{1,i}=uiEnterName(fmriprep_table_names{i,1},'Enter featureName');
    end
end

%% set parcellations to load
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
    ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
    for i = 1:length(ParcelNames)
        ParcelNames{i,2}=uiEnterName(ParcelNames{i,1},'Enter parcel abbreviation');
    end     
end

%%Compile Behavioral Data
if isempty(TaskVarNames)
    SingleSelect=1; %Allows only a single value to be selected.
    if isempty(IncludeTaskVars)
        [IncludeTaskVars] = uiNameSelect({'Yes','No'},'Include task vars as covariates?',SingleSelect);
        if strcmpi(IncludeTaskVars,'Yes')
            IncludeTaskVars=1;
        else
            IncludeTaskVars=0;
        end
    end
elseif strcmpi(TaskVarNames,'NA_EMPTY')
    IncludeTaskVars=0;
    TaskVarNames=[];
else
    IncludeTaskVars=1;
end

AllTaskVarNames=[];

for i = 1:NumGroups
    try
    [tempBehavior_data_all] = CompileND(ExperimentsDirs{i,1},fmriprep_tables{i,1},...
        'LoadVarName','beh_vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','beh_vars',...
        'AnalysisName','Overall',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    catch
    [tempBehavior_data_all] = CompileND(ExperimentsDirs{i,1},fmriprep_tables{i,1},...
        'LoadVarName','SS_Vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','SSVars',...
        'AnalysisName','SSVars_All',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    end
    AllTaskVarNames=unique([AllTaskVarNames;tempBehavior_data_all.Properties.VariableNames(:)]);
end
if isempty(TaskVarNames)
    TaskVarNames=uiNameSelect([{'none'};AllTaskVarNames],'Select behavioral variables',0);
end
if ~isempty(TaskVarNames)
    if ~iscell(beh_var_names)
        TaskVarNames={TaskVarNames};
    end
    if any(strcmpi(TaskVarNames,'none')) || any(strcmpi(TaskVarNames,'NA_EMPTY'))
        TaskVarNames=[];
    end
else
    TaskVarNames=[];
end
    
AnalysisParameters.TaskVarNames=TaskVarNames;


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
if isempty(SSVarNames)
    if IncludeSSVars==1
        AllSSVarNames=[];
        for i = 1:NumGroups
            [tempBehavior_data_all] = CompileND(ExperimentsDirs{i,1},fmriprep_tables{i,1},...
                'LoadVarName','SS_Vars',...
                'LoadVarFormat','Table',...
                'DataType','Other',...
                'AnalysisType','SSVars',...
                'AnalysisName','SSVars_All',...
                'SubjectOrRun',SubjectOrRun,...
                'TableOrArray','Table',...
                'ParcelName','None');
            AllSSVarNames=unique([AllSSVarNames;tempBehavior_data_all.Properties.VariableNames(:)]);
        end
        if isempty(SSVarNames)
            SSVarNames=uiNameSelect([{'none'};AllSSVarsNames],'Select behavioral variables',0);
        end
    else
        SSVarNames=[];
    end
end
for i = 1:NumGroups
    [AllTables{i,1},BehTable{i,1},ssFeatureTable{i,1},selectInd{i,1},use_fmriprep_table{i,1}] = compileMultiFeatureTable(fmriprep_tables{i,1},ExperimentsDirs{i,1},ParcelNames,featureNames,...
        'IncludeTaskVars',IncludeTaskVars,'IncludeSSVars',IncludeSSVars,'SSVarNames',SSVarNames,'TaskVarNames',TaskVarNames,'subIDNames',subIDNames{i,1},'SubjectOrRun',SubjectOrRun);
end

%% Define LME model
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

end

