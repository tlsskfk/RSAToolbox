function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeSummaryStats_IndvDiff(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Overwrite previously saved files (default is no or 0; yes = 1)
[ExperimentsDir] = VariableSetter('ExperimentsDir',[],varargin);
%Overwrite previously saved files (default is no or 0; yes = 1)
[fmriprep_table] = VariableSetter('fmriprep_table',[],varargin);
%Overwrite previously saved files (default is no or 0; yes = 1)
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[LoadVarName] = VariableSetter('LoadVarName',[],varargin);
[LoadVarFormat] = VariableSetter('LoadVarFormat',[],varargin);
[DataType] = VariableSetter('DataType',[],varargin);
[CompileType] = VariableSetter('CompileType',[],varargin); %TableOrArray
[CompileFormat] = VariableSetter('CompileFormat',[],varargin); %VertOrMat
[ndim] = VariableSetter('ndim',[],varargin);

AnalysisType = VariableSetter('AnalysisName',[],varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

All_Analyses=cell(1);

[AnalysisParameters.filePaths,AnalysisParameters.AnalysisType,AnalysisParameters.AnalysisName,AnalysisParameters.ExperimentsDir,AnalysisParameters.fmriprep_table,AnalysisParameters.fmriprep_table_name,AnalysisParameters.SubjectOrRun] = BIDsDirSearch(ExperimentsDir,fmriprep_table,...
    'SubjectOrRun',SubjectOrRun,...
    'AnalysisType',AnalysisType,...
    'AnalysisName',AnalysisName,...
    'TitleTextName',['Select files to compile']);
ParcelNames=AnalysisParameters.filePaths.Properties.VariableNames;
AnalysisParameters.ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run',0);
AnalysisParameters.ParcelNames=AnalysisParameters.ParcelNames(:);
[AnalysisParameters.LoadVarName,AnalysisParameters.DataType,AnalysisParameters.LoadVarFormat,AnalysisParameters.CompileType,AnalysisParameters.CompileFormat,AnalysisParameters.ndim] = GetFileVariableNames(filePaths,'LoadVarName',LoadVarName,...
    'LoadVarFormat',LoadVarFormat,...
    'DataType',DataType,...
    'CompileType',CompileType,...
    'CompileFormat',CompileFormat,...
    'ndim',ndim);

SetAnalysesLoop=0;
count=1;
while SetAnalysesLoop == 0
    All_Analyses{count,1}=AnalysisParameters;
    [AddAnalysis] = uiNameSelect({'Yes','No'},'Run additional analysis?',SingleSelect);   
    if strcmpi(AddAnalysis,'Yes')    
        [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect); 
    else
        SetAnalysesLoop=1;
        continue
    end
    for i = 1:length(EditParams)
        AnalysisParameters.(EditParams{i,1})=[];
    end   
    [AnalysisParameters.filePaths,AnalysisParameters.AnalysisType,AnalysisParameters.AnalysisName,AnalysisParameters.ExperimentsDir,AnalysisParameters.fmriprep_table,AnalysisParameters.fmriprep_table_name,AnalysisParameters.SubjectOrRun] = BIDsDirSearch(AnalysisParameters.ExperimentsDir,AnalysisParameters.fmriprep_table,...
        'SubjectOrRun',AnalysisParameters.SubjectOrRun,...
        'AnalysisType',AnalysisParameters.AnalysisType,...
        'AnalysisName',AnalysisParameters.AnalysisName,...
        'TitleTextName',['Select files to compile']);
    if isempty(AnalysisParameters.ParcelNames)
        ParcelNames=filePaths.Properties.VariableNames;
        AnalysisParameters.ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run',0);
        AnalysisParameters.ParcelNames=AnalysisParameters.ParcelNames(:);
    end
    [AnalysisParameters.LoadVarName,AnalysisParameters.DataType,AnalysisParameters.LoadVarFormat,AnalysisParameters.CompileType,AnalysisParameters.CompileFormat,AnalysisParameters.ndim] = GetFileVariableNames(AnalysisParameters.filePaths,'LoadVarName',AnalysisParameters.LoadVarName,...
        'LoadVarFormat',AnalysisParameters.LoadVarFormat,...
        'DataType',AnalysisParameters.DataType,...
        'CompileType',AnalysisParameters.CompileType,...
        'CompileFormat',AnalysisParameters.CompileFormat,...
        'ndim',AnalysisParameters.ndim);    
    count=count+1;
end


for i = 1:size(All_Analyses,1)
    %%Define Variables
    AnalysisParameters=All_Analyses{i,1};
    ExperimentsDir=AnalysisParameters.ExperimentsDir;
    fmriprep_table=AnalysisParameters.fmriprep_table;
    fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    AnalysisType=AnalysisParameters.AnalysisType;
    AnalysisName=AnalysisParameters.AnalysisName;
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    ParcelNames=AnalysisParameters.ParcelNames;
    LoadVarName=AnalysisParameters.LoadVarName;
    DataType=AnalysisParameters.DataType;
    LoadVarFormat=AnalysisParameters.LoadVarFormat;
    CompileType=AnalysisParameters.CompileType;
    CompileFormat=AnalysisParameters.CompileFormat;
    ndim=AnalysisParameters.ndim;
    filePaths=AnalysisParameters.filePaths;
    %%Specify group results directory and create if it doesn't exist
    ExperimentsDir=strrep(ExperimentsDir,'\','/');
    GroupDir=strrep(ExperimentsDir,'/Experiments/','');
    GroupDir=[GroupDir,'/GroupAnalysis/'];
    GroupDir=strrep(GroupDir,'//','/');
    FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
    BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
    for j = 1:size(ParcelNames,1)
        ParcelName=ParcelNames{j,1};
        [compiledData,compiled_fmriprep_table,selectind,DataLabels] = CompileND(ExperimentsDir,...
            AnalysisParameters.fmriprep_table,...
            'ParcelName',ParcelName,...
            'AnalysisType',AnalysisType,...
            'AnalysisName',AnalysisName,...
            'LoadVarFormat',LoadVarFormat,...
            'LoadVarName',LoadVarName,...
            'DataType',DataType,...
            'SubjectOrRun',SubjectOrRun,...
            'TableOrArray',CompileType,...
            'VertOrMat',CompileFormat);
        numNs=height(compiled_fmriprep_table);
        if ndim==3
            if size(compiledData,1)~=size(compiledData,2)
                disp(['unknown format! skipping!'])
                continue
            end
            compiledData=mat2uppertriuvectormat(compiledData)';
        end
        if size(compiledData,1) == size(compiledData,2)
            disp(['warning! check data orientation']);
        elseif size(compiledData,2) == numNs
            compiledData=compiledData';
        end
        
    end
    
end



end

