function [AllTables,BehTable,ssFeatureTable,selectInd,use_fmriprep_table] = compileMultiFeatureTable(fmriprep_table,ExperimentsDir,ParcelNames,featureNames,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end
[IncludeTaskVars] = VariableSetter('IncludeTaskVars',[],varargin); %1 = yes 
[TaskVarNames] = VariableSetter('TaskVarNames',[],varargin);
[IncludeSSVars] = VariableSetter('IncludeSSVars',[],varargin);
[SSVarNames] = VariableSetter('SSVarNames',[],varargin);
[subIDNames] = VariableSetter('subIDNames',[],varargin);
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
ssFeatureTable=[];
AllTables=cell(size(ParcelNames,1),length(featureNames));
IndDiffTables=cell(size(ParcelNames,1),length(featureNames));
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
if isempty(subIDNames)
    SingleSelect=1; %Allows only a single value to be selected.
    [subIDNames] = uiNameSelect(fmriprep_table.Properties.VariableNames,'Select variable for SS names:',SingleSelect);
end
if iscell(subIDNames)
    subIDNames=subIDNames{1,1};
end
subIDNames=fmriprep_table.(subIDNames);
if strcmpi(SubjectOrRun,'Run')
    try
        subIDNames=strcat(subIDNames,num2str(fmriprep_table.session),num2str(fmriprep_table.run));
    catch
        subIDNames=strcat(subIDNames,num2str(fmriprep_table.run));
    end
end
%%Compile Task variables
SingleSelect=1; %Allows only a single value to be selected.
if isempty(IncludeTaskVars)
    [IncludeTaskVars] = uiNameSelect({'Yes','No'},'Include task vars?',SingleSelect);
    if strcmpi(IncludeTaskVars,'Yes')
        IncludeTaskVars=1;
    else
        IncludeTaskVars=0;
    end
end
selectInd=ones(height(fmriprep_table),1);
if IncludeTaskVars==1
    [TaskVars_data_all,~,tempInd] = CompileND(ExperimentsDir,fmriprep_table,...
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
    TaskVars_data_all.Properties.RowNames=subIDNames(tempInd==1);
    if isempty(ssFeatureTable)
        ssFeatureTable=TaskVars_data_all;
    else
        ssFeatureTable=join(ssFeatureTable,TaskVars_data_all,'keys','RowNames');
    end
    selectInd=selectInd.*tempInd;
    clearvars TaskVars_data_all
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
    [SSVars_data_all,~,tempInd] = CompileND(ExperimentsDir,fmriprep_table,...
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
    if isempty(ssFeatureTable)
        ssFeatureTable=SSVars_data_all;
    else
        try
            ssFeatureTable=join(ssFeatureTable,SSVars_data_all,'keys','RowNames');
        catch
            ssFeatureTable=join(SSVars_data_all,ssFeatureTable,'keys','RowNames');
        end
    end
    selectInd=selectInd.*tempInd;
    clearvars SSVars_data_all
end
BehTable=ssFeatureTable;
numParcels=size(ParcelNames,1);
numFeatureNames=length(featureNames);
if size(ParcelNames,2)==2
    ParcelAbbrvs=ParcelNames(:,2);
else
    ParcelAbbrvs=ParcelNames(:,1);
end

for FeatNum = 1:numFeatureNames
    featureName=featureNames{FeatNum};
    [RunVarTypes] = uiNameSelect({'Nodes','Edges','Both'},['Select feature type for ',featureName],SingleSelect);
    if strcmpi(RunVarTypes,'Nodes') || strcmpi(RunVarTypes,'Both')
        for parcelNum=1:numParcels
%             try
                parcelName=ParcelNames{parcelNum,1};
                disp(['Loading ',featureName,'-Nodes for parcellation: ',parcelName,'.']);
                parcelAbbrv=ParcelAbbrvs{parcelNum,1};
                if parcelNum == 1
                    CompileLoadParams.LoadVarName=[];
                    CompileLoadParams.LoadVarFormat='Table';
                    CompileLoadParams.DataType='ByParcellation';
                    CompileLoadParams.AnalysisType=[];
                    CompileLoadParams.AnalysisName=[];
                    CompileLoadParams.SubjectOrRun=SubjectOrRun;
                    CompileLoadParams.TableOrArray='Table';
                    CompileLoadParams.RowCompile=[];
                end 
                [tempDataTable,~,tempInd,DataLabels,~,~,~,CompileLoadParams] = CompileND(ExperimentsDir,fmriprep_table,...
                    'LoadVarName',CompileLoadParams.LoadVarName,...
                    'LoadVarFormat',CompileLoadParams.LoadVarFormat,...
                    'DataType',CompileLoadParams.DataType,...
                    'AnalysisType',CompileLoadParams.AnalysisType,...
                    'AnalysisName',CompileLoadParams.AnalysisName,...
                    'SubjectOrRun',CompileLoadParams.SubjectOrRun,...
                    'TableOrArray',CompileLoadParams.TableOrArray,...
                    'ParcelName',parcelName,...
                    'RowCompile',CompileLoadParams.RowCompile);           
                if ~istable(tempDataTable)
                    if ndims(tempDataTable)==3
                        tempDataTable=mat2uppertriuvectormat(tempDataTable)';
                        DataLabels=labels2uppertriuvectorlabels(DataLabels);
                    end
                    tempDataTable=array2table(tempDataTable,'RowNames',subIDNames(tempInd==1),'VariableNames',DataLabels);
                elseif isempty(tempDataTable.Properties.RowNames)
                    try
                        tempDataTable.Properties.RowNames=subIDNames(tempInd==1);
                    end
                end   
                tempDataTable.Properties.VariableNames=join([repmat({featureName},[length(tempDataTable.Properties.VariableNames),1]),tempDataTable.Properties.VariableNames(:),repmat({parcelAbbrv},[length(tempDataTable.Properties.VariableNames),1])],'_');
                if isempty(ssFeatureTable)
                    ssFeatureTable=tempDataTable;
                else
                    tempDataTable=tempDataTable(ismember(tempDataTable.Properties.RowNames,ssFeatureTable.Properties.RowNames),:);
                    try
                        ssFeatureTable=join(ssFeatureTable,tempDataTable,'keys','RowNames');
                    catch
                        ssFeatureTable=join(tempDataTable,ssFeatureTable,'keys','RowNames');
                    end
                end 
                selectInd=selectInd.*tempInd;
                AllTables{parcelNum,FeatNum}=tempDataTable;
%             catch
%                 disp('No data--skipping')
%             end
        end        
    end        
    if strcmpi(RunVarTypes,'Edges') || strcmpi(RunVarTypes,'Both')
        for parcelNum=1:numParcels
            parcelName=ParcelNames{parcelNum,1};
            disp(['Loading ',featureName,'-Edges for parcellation: ',parcelName,'.']);
            parcelAbbrv=ParcelAbbrvs{parcelNum,1};
            if parcelNum == 1
                CompileLoadParams.LoadVarName=[];
                CompileLoadParams.LoadVarFormat='Table';
                CompileLoadParams.DataType='ByParcellation';
                CompileLoadParams.AnalysisType=[];
                CompileLoadParams.AnalysisName=[];
                CompileLoadParams.SubjectOrRun=SubjectOrRun;
                CompileLoadParams.TableOrArray='Table';
                CompileLoadParams.RowCompile=[];
            end 
            [tempDataTable,~,tempInd,DataLabels,~,~,~,CompileLoadParams] = CompileND(ExperimentsDir,fmriprep_table,...
                'LoadVarName',CompileLoadParams.LoadVarName,...
                'LoadVarFormat',CompileLoadParams.LoadVarFormat,...
                'DataType',CompileLoadParams.DataType,...
                'AnalysisType',CompileLoadParams.AnalysisType,...
                'AnalysisName',CompileLoadParams.AnalysisName,...
                'SubjectOrRun',CompileLoadParams.SubjectOrRun,...
                'TableOrArray',CompileLoadParams.TableOrArray,...
                'ParcelName',parcelName,...
                'RowCompile',CompileLoadParams.RowCompile);  
            if ~istable(tempDataTable)
                if ndims(tempDataTable)==3
                    tempDataTable=mat2uppertriuvectormat(tempDataTable)';
                    DataLabels=labels2uppertriuvectorlabels(DataLabels);
                end
                DataLabels=truncateStrCell(DataLabels,namelengthmax);
                tempDataTable=array2table(tempDataTable,'RowNames',subIDNames(tempInd==1),'VariableNames',DataLabels);
            elseif isempty(tempDataTable.Properties.RowNames)
                try
                    tempDataTable.Properties.RowNames=subIDNames(tempInd==1);
                end
            end   
            tempDataTable.Properties.VariableNames=truncateStrCell(join([repmat({featureName},[length(tempDataTable.Properties.VariableNames),1]),tempDataTable.Properties.VariableNames(:),repmat({parcelAbbrv},[length(tempDataTable.Properties.VariableNames),1])],'_'),namelengthmax);
            
            if isempty(ssFeatureTable)
                ssFeatureTable=tempDataTable;
            else
                tempDataTable=tempDataTable(ismember(tempDataTable.Properties.RowNames,ssFeatureTable.Properties.RowNames),:);
                try
                    ssFeatureTable=join(ssFeatureTable,tempDataTable,'keys','RowNames');
                catch
                    ssFeatureTable=join(tempDataTable,ssFeatureTable,'keys','RowNames');
                end
            end 
            selectInd=selectInd.*tempInd;
            AllTables{parcelNum,FeatNum}=tempDataTable;
        end               
    end
end
use_fmriprep_table=fmriprep_table(selectInd==1,:);
useIDs = subIDNames(selectInd==1,:);
if ~isempty(BehTable)
    BehTable=BehTable(useIDs,:);
end
for FeatNum = 1:numFeatureNames
    for parcelNum=1:numParcels
        AllTables{parcelNum,FeatNum}=AllTables{parcelNum,FeatNum}(useIDs,:);
    end
end


end

