function [compiledData,compiled_fmriprep_table,selectind,ROILabels] = CompileRSMs(ExperimentsDir,fmriprep_table,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[AnalysisType] = VariableSetter('AnalysisType',['RSMs'],varargin);
[LoadVarName] = VariableSetter('LoadVarName',['RSMs'],varargin);

[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[ParcelName] = VariableSetter('ParcelName',[],varargin);
[VertOrMat] = VariableSetter('VertOrMat',[],varargin);
[TableOrArray] = VariableSetter('TableOrArray',[],varargin);


if isempty(VertOrMat)
    SingleSelect=1; %Allows only a single value to be selected.
    [VertOrMat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',SingleSelect);
end
if strcmpi(VertOrMat,'Matrix')
    TableOrArray='Array';
end

if isempty(TableOrArray)
    SingleSelect=1; %Allows only a single value to be selected.
    [TableOrArray] = uiNameSelect({'Table','Array'},'Table or array form?',SingleSelect);
end

[filePaths,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',AnalysisType,'AnalysisName',AnalysisName,'SubjectOrRun',SubjectOrRun);
if isempty(parcelName)
    ParcelNames=filePaths.Properties.VariableNames;
    [ParcelName] = uiNameSelect(ParcelNames,'Select parcel to compile:' ,1);
end
filePaths=filePaths.(ParcelName);
compiledData=[];
selectind=zeros(height(fmriprep_table),1);
getLabels=0;
for i =1:length(filePaths)
    try
        LoadData=load(filePaths{i,1},LoadVarName);
        LoadData=LoadData.(LoadVarName);
        if getLabels==0            
            if strcmpi(VertOrMat,'Vertical')
                ROILabels.ROILabels=LoadData.Properties.VariableNames;
                [ROILabels.labelPairs1,ROILabels.labelPairs2,ROILabels.labelPairsCell] = labels2uppertriuvectorlabels( ROILabels.ROILabels );
            else
                ROILabels=LoadData.Properties.VariableNames;
            end
            getLabels=1;
        end
        LoadData=table2array(LoadData);
        if strcmpi(VertOrMat,'Vertical')
            LoadData=mat2uppertriuvectormat(LoadData)';
            if strcmpi(TableOrArray,'Table')
                LoadData=array2table(LoadData,'VariableNames',ROILabels.labelPairs2);
            end
            compiledData=cat(2,compiledData,LoadData);
        else
            compiledData=cat(3,compiledData,LoadData);
        end
    catch
        continue
    end
    selectind(i,1)=1;
    
end

if strcmpi(TableOrArray,'Table') 
    subNames=fmriprep_table.sub(selectind==1,:);
    compiledData=[subNames,compiledData];
end
compiled_fmriprep_table=fmriprep_table(selectind,:);

end

