function [compiledData,compiled_fmriprep_table,selectind,DataLabels,ExperimentsDir,fmriprep_table,fmriprep_table_name,LoadParams] = CompileND(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% [compiledData,compiled_fmriprep_table,selectind,DataLabels,ExperimentsDir,fmriprep_table,fmriprep_table_name,LoadParams] = CompileND(ExperimentsDir,fmriprep_table,...
%     'AnalysisType',AnalysisType,...
%     'LoadVarName',LoadVarName,...
%     'LoadVarFormat',LoadVarFormat,...
%     'DataType',DataType,...
%     'AnalysisName',AnalysisName,...
%     'SubjectOrRun',SubjectOrRun,...
%     'ParcelName',ParcelName,...
%     'TableOrArray',TableOrArray,...
%     'VertOrMat','Matrix',...
%     'RowCompile',RowCompile,...
%     'ColumnCompile','All');
    
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end

[AnalysisType] = VariableSetter('AnalysisType',[],varargin);
[LoadVarName] = VariableSetter('LoadVarName',[],varargin);

[LoadVarFormat] = VariableSetter('LoadVarFormat',[],varargin);
[DataType] = VariableSetter('DataType',[],varargin);
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[ParcelName] = VariableSetter('ParcelName',[],varargin);
[TableOrArray] = VariableSetter('TableOrArray',[],varargin);
[VertOrMat] = VariableSetter('VertOrMat',['Matrix'],varargin);
[RowCompile] = VariableSetter('RowCompile',[],varargin);
[ColumnCompile] = VariableSetter('ColumnCompile',['All'],varargin);

if isempty(VertOrMat)
    SingleSelect=1; %Allows only a single value to be selected.
    [VertOrMat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',SingleSelect);
end
% if strcmpi(VertOrMat,'Matrix')
%     TableOrArray='Array';
% end
if isempty(DataType)
    SingleSelect=1; %Allows only a single value to be selected.
    [DataType] = uiNameSelect({'ByParcellation','ByCoordinate','Other'},'How is the data organized?',SingleSelect);
end
if ~strcmpi(DataType,'ByParcellation')
    ParcelName='None';
end

if strcmpi(LoadVarFormat,'Table')
    DataLabels=[];
else
    DataLabels=['None'];
end

if isempty(TableOrArray)
    SingleSelect=1; %Allows only a single value to be selected.
    [TableOrArray] = uiNameSelect({'Table','Array'},'Compile in table or array form?',SingleSelect);
end

[filePaths,AnalysisType,AnalysisName,~,~,~,SubjectOrRun] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',AnalysisType,'AnalysisName',AnalysisName,'SubjectOrRun',SubjectOrRun);

if isempty(ParcelName)
    ParcelNames=filePaths.Properties.VariableNames;
    [ParcelName] = uiNameSelect(ParcelNames,'Select parcel to compile:' ,1);
end
if strcmpi(ParcelName,'None')
    if length(AnalysisName) > namelengthmax
        filePaths=filePaths.(AnalysisName(1,1:namelengthmax));
    else
        filePaths=filePaths.(AnalysisName);
    end
else
    try
        filePaths=filePaths.(ParcelName);
    catch
        filePaths=filePaths.(ParcelName{1,1});
    end
    if strcmpi(DataLabels,'None')
        try
            load(['Parcellations/',ParcelName],'UseLabels')
            if strcmpi(VertOrMat,'Vertical')
                DataLabels.ROILabels=UseLabels;
                [DataLabels.labelPairs1,DataLabels.labelPairs2,DataLabels.labelPairsCell] = labels2uppertriuvectorlabels( DataLabels.ROILabels );   
            else
                DataLabels=UseLabels(:);
            end
        catch
            DataLabels=[];
        end
    end
end
if isempty(LoadVarName)
    for i = 1:length(filePaths)
        try
            variableInfo = who('-file', filePaths{i,1});
            LoadVarName=unique([LoadVarName;variableInfo]);
        catch
            continue
        end
    end
    [LoadVarName] = uiNameSelect(LoadVarName,'Select variable to load:' ,1);
end
if iscell(LoadVarName)
    LoadVarName=LoadVarName{1,1};
end
compileDim=[];
compiledData=[];
selectind=zeros(height(fmriprep_table),1);
for i =1:length(filePaths)
    try
        LoadData=load(filePaths{i,1},LoadVarName);
        LoadData=LoadData.(LoadVarName);
        if isempty(LoadVarFormat)
            if istable(LoadData)
                LoadVarFormat='Table';                
            elseif iscell(LoadData)
                LoadVarFormat='Cell';
            elseif isnumeric(LoadData)
                LoadVarFormat='Array';
            else
                disp('Error: unknown data format')
                break
            end
            if strcmpi(LoadVarFormat,'Table')
                DataLabels=[];
            end
        end
        %if strcmpi(LoadVarFormat,'Table') && strcmpi(TableOrArray,'Table')
        if strcmpi(LoadVarFormat,'Table')
            if sum(ismember(LoadData.Properties.RowNames,LoadData.Properties.VariableNames)) ~= length(LoadData.Properties.VariableNames)
                if height(LoadData)>1     
                    if isempty(RowCompile) && ~isempty(LoadData.Properties.RowNames)
                        RowNames=LoadData.Properties.RowNames;
                        RowCompile=uiNameSelect([{'All'};RowNames(:)],'Select row variable:',1);
                    end
                    if ~isempty(RowCompile) && ~strcmpi(RowCompile,'All')
                        LoadData=LoadData(RowCompile,:);
                        LoadData.Properties.RowNames={};
                    end
                end
                if size(LoadData,2) > 1
                    if isempty(ColumnCompile) && ~isempty(LoadData.Properties.VariableNames)
                        ColumnNames=LoadData.Properties.VariableNames;
                        ColumnCompile=uiNameSelect([{'All'};ColumnNames(:)],'Select column variable:',1);
                    end
                    if ~isempty(ColumnCompile) && ~strcmpi(ColumnCompile,'All')
                        LoadData=LoadData(ColumnNames,:);
                        LoadData.Properties.VariableNames={};
                    end  
                end
            end
        end
        if isempty(DataLabels) && strcmpi(LoadVarFormat,'Table')            
            if strcmpi(VertOrMat,'Vertical')
                DataLabels.ROILabels=LoadData.Properties.VariableNames;
                [DataLabels.labelPairs1,DataLabels.labelPairs2,DataLabels.labelPairsCell] = labels2uppertriuvectorlabels( DataLabels.ROILabels );
            else
                DataLabels=LoadData.Properties.VariableNames;
            end
        end

        if strcmpi(LoadVarFormat,'Table') && strcmpi(VertOrMat,'Vertical')
            try
                LoadData=table2array(LoadData);
            catch
                TableOrArray='Table';
                VertOrMat='Matrix';               
            end
        end
        if strcmpi(VertOrMat,'Vertical')
            LoadData=mat2uppertriuvectormat(LoadData)';
            if strcmpi(TableOrArray,'Table')
                LoadData=array2table(LoadData,'VariableNames',DataLabels.labelPairs2);
            end
            compiledData=cat(2,compiledData,LoadData);
        else
            if strcmpi(TableOrArray,'Array') && istable(LoadData)
                LoadData=table2array(LoadData);
            end
            if isempty(compileDim)
                if isvector(LoadData)
                    if size(LoadData,1)==1
                        compileDim=1;
                    else
                        compileDim=2;
                    end                       
                else
                    compileDim=ndims(LoadData)+1;
    
                end         
            end
            if compileDim > 2 && istable(LoadData)
                LoadData=table2array(LoadData); 
                TableOrArray='Array';
            end
            compiledData=cat(compileDim,compiledData,LoadData);
        end
    catch
        continue
    end
    selectind(i,1)=1;
    
end

if strcmpi(TableOrArray,'Table')
    try
        subNames=fmriprep_table.orig_IDs(selectind==1,:);
    catch    
        subNames=fmriprep_table.sub(selectind==1,:);
    end
    try
        compiledData.Properties.RowNames=subNames;
    end
end

compiled_fmriprep_table=fmriprep_table(selectind==1,:);
LoadParams.AnalysisType=AnalysisType;
LoadParams.LoadVarName=LoadVarName;
LoadParams.LoadVarFormat=LoadVarFormat;
LoadParams.DataType=DataType;
LoadParams.AnalysisName=AnalysisName;
LoadParams.SubjectOrRun=SubjectOrRun;
LoadParams.ParcelName=ParcelName;
LoadParams.TableOrArray=TableOrArray;
LoadParams.VertOrMat=VertOrMat;
LoadParams.ColumnCompile=ColumnCompile;
LoadParams.RowCompile=RowCompile;
end

