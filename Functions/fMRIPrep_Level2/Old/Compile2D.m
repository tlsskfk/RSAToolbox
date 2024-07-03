function [compiledData,compiled_fmriprep_table,selectind,DataLabels,ExperimentsDir,fmriprep_table,fmriprep_table_name] = Compile2D(ExperimentsDir,fmriprep_table,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
[VertOrMat] = VariableSetter('VertOrMat',[],varargin);

if isempty(VertOrMat)
    SingleSelect=1; %Allows only a single value to be selected.
    [VertOrMat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',SingleSelect);
end
if strcmpi(VertOrMat,'Matrix')
    TableOrArray='Array';
end
if isempty(DataType)
    SingleSelect=1; %Allows only a single value to be selected.
    [DataType] = uiNameSelect({'ByParcellation','ByCoordinate','Other'},'How is the data organized?',SingleSelect);
end
if ~strcmpi(DataType,'ByParcellation')
    ParcelName='None';
end
if isempty(LoadVarFormat)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadVarFormat] = uiNameSelect({'Table','Array','Cell'},'What is the format of the data to compile?',SingleSelect);
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

[filePaths,~,AnalysisName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',AnalysisType,'AnalysisName',AnalysisName,'SubjectOrRun',SubjectOrRun);

if isempty(ParcelName)
    ParcelNames=filePaths.Properties.VariableNames;
    [ParcelName] = uiNameSelect(ParcelNames,'Select parcel to compile:' ,1);
end
if strcmpi(ParcelName,'None')
    filePaths=filePaths.(AnalysisName);
else
    filePaths=filePaths.(ParcelName);
    if strcmpi(DataLabels,'None')
        load(['Parcellations/',ParcelName],'UseLabels')
        if strcmpi(VertOrMat,'Vertical')
            DataLabels.ROILabels=UseLabels;
            [DataLabels.labelPairs1,DataLabels.labelPairs2,DataLabels.labelPairsCell] = labels2uppertriuvectorlabels( DataLabels.ROILabels );   
        else
            DataLabels=UseLabels(:);
        end
    end
end

compiledData=[];
selectind=zeros(height(fmriprep_table),1);
for i =1:length(filePaths)
    try
        LoadData=load(filePaths{i,1},LoadVarName);
        LoadData=LoadData.(LoadVarName);
        if isempty(DataLabels)            
            if strcmpi(VertOrMat,'Vertical')
                DataLabels.ROILabels=LoadData.Properties.VariableNames;
                [DataLabels.labelPairs1,DataLabels.labelPairs2,DataLabels.labelPairsCell] = labels2uppertriuvectorlabels( DataLabels.ROILabels );
            else
                DataLabels=LoadData.Properties.VariableNames;
            end
        end
        if strcmpi(LoadVarFormat,'Table')
            LoadData=table2array(LoadData);
        end
        if strcmpi(VertOrMat,'Vertical')
            LoadData=mat2uppertriuvectormat(LoadData)';
            if strcmpi(TableOrArray,'Table')
                LoadData=array2table(LoadData,'VariableNames',DataLabels.labelPairs2);
            end
            compiledData=cat(1,compiledData,LoadData);
        else
            compiledData=cat(3,compiledData,LoadData);
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
%    compiledData.Properties.RowNames=subNames;
end
compiled_fmriprep_table=fmriprep_table(selectind==1,:);

end

