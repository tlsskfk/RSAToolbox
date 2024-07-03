function [compiledData,compiled_fmriprep_table,selectind,ROILabels] = CompileRF(ExperimentsDir,fmriprep_table,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[AnalysisType] = VariableSetter('AnalysisType',['RCA'],varargin);
[LoadVarName] = VariableSetter('LoadVarName',['RCA_RF'],varargin);

[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[ParcelName] = VariableSetter('ParcelName',[],varargin);
[TableOrArray] = VariableSetter('TableOrArray',[],varargin);

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
            ROILabels=LoadData.Properties.VariableNames;
            getLabels=1;
        end
        if strcmpi(TableOrArray,'Array')
            compiledData=table2array(compiledData);
        end
        compiledData=cat(2,compiledData,LoadData);
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

