function [compiledData,DataInfo] = uiDataCompiler(varargin)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[ExperimentsDir] = VariableSetter('ExperimentsDir',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
[fmriprep_table] = VariableSetter('fmriprep_table',[],varargin);
% [AnalysisName] = VariableSetter('AnalysisName',[],varargin);
% [AnalysisType] = VariableSetter('AnalysisType',[],varargin);
global uiDCParameters;
UIFig=uifigure('Position',[500,50,500,500],'Visible','On');
ListBox=uilistbox(UIFig,'Position',[20,60,260,250],'Items',{},'Multiselect','off');
if isempty(ExperimentsDir)
    ExpDir=uilabel(UIFig,'Position',[130,500-50,380,50],'Text','Press button to select experiment dir','HorizontalAlignment','left');  
    SelectExpDirButton=uibutton(UIFig,'Position',[20,500-45,100,40],'Visible','On','Text','SelectDir','ButtonPushedFcn',@(SelectExpDirButton,event) GetExpDirPush(uiDCParameters,ExpDir) );
else
    ExpDir=uilabel(UIFig,'Position',[130,500-50,380,50],'Text',ExperimentsDir,'HorizontalAlignment','left');  
    SelectExpDirButton=uibutton(UIFig,'Position',[20,500-45,100,40],'Visible','Off','Text','SelectDir','ButtonPushedFcn',@(SelectExpDirButton,event) GetExpDirPush(uiDCParameters,ExpDir) );
    uiDCParameters.ExperimentsDir=ExperimentsDir;
end

if ~isempty(fmriprep_table)
    uiDCParameters.fmriprep_table=fmriprep_table;
else
    if isempty(fmriprep_table_name)
        FmriPrepTableName=uilabel(UIFig,'Position',[130,500-100,380,50],'Text','Press button to select fmriprep table','HorizontalAlignment','left');  
        SelectFmriPrepTableName=uibutton(UIFig,'Position',[20,500-95,100,40],'Visible','On','Text','Select table','ButtonPushedFcn',@(SelectFmriPrepTableName,event) ...
            GetTableNamePush(uiDCParameters,FmriPrepTableName) );
    else
        FmriPrepTableName=uilabel(UIFig,'Position',[130,500-100,380,50],'Text',fmriprep_table_name,'HorizontalAlignment','left');  
        SelectFmriPrepTableName=uibutton(UIFig,'Position',[20,500-95,100,40],'Visible','Off','Text','Select table','ButtonPushedFcn',@(SelectFmriPrepTableName,event) ...
            GetTableNamePush(uiDCParameters,FmriPrepTableName) );
        uiDCParameters.fmriprep_table_name=fmriprep_table_name;
        load(['fmriprep_table/',fmriprep_table_name],'fmriprep_table');
        uiDCParameters.fmriprep_table=fmriprep_table;
    end
end
if ~isempty(fmriprep_table_name)
    FmriPrepTableName=uilabel(UIFig,'Position',[130,500-100,380,50],'Text',fmriprep_table_name,'HorizontalAlignment','left');  
end

SubjectOrRun = uiswitch(UIFig,'Position',[40,500-135,380,30],'Items',{'Run','SS'},'Value','Run');
AnalysisTypeLabel=uilabel(UIFig,'Position',[320,500-165,180,30],'Text','Analysis Type: ','HorizontalAlignment','left');  
AnalysisNameLabel=uilabel(UIFig,'Position',[320,500-200,180,30],'Text','Analysis Name: ','HorizontalAlignment','left');  
ParcellationLabel=uilabel(UIFig,'Position',[320,500-235,180,30],'Text','Parcellation: ','HorizontalAlignment','left');  
VariableLabel=uilabel(UIFig,'Position',[320,500-270,180,30],'Text','Load Variable: ','HorizontalAlignment','left');  

DataType = uiswitch(UIFig,'Position',[400,500-310,180,30],'Visible','Off','Items',{'ByParcellation','Other'},'Value','Other');
TableOrArray = uiswitch(UIFig,'Position',[400,500-350,180,30],'Visible','Off','Items',{'Table','Array'},'Value','Array');
VertOrMat = uiswitch(UIFig,'Position',[400,500-390,180,30],'Visible','Off','Items',{'Vertical','Matrix'},'Value','Matrix');




CompileDataButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','Off','Text','Compile Data','ButtonPushedFcn',@(CompileDataButton,event) CompileDataPush(uiDCParameters,UIFig,VertOrMat,TableOrArray,DataType) );

SelectVarButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','Off','Text','Select Variable','ButtonPushedFcn',@(SelectVarButton,event) ...
    CompileVarInfo(uiDCParameters,ListBox,VariableLabel,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat) );
SelectParcelButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','Off','Text','Select Parcel','ButtonPushedFcn',@(SelectParcelButton,event) ...
    CompileVariables(uiDCParameters,ListBox,ParcellationLabel,VariableLabel,SelectParcelButton,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat) );
SelectNameButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','Off','Text','Select Name','ButtonPushedFcn',@(SelectNameButton,event) ...
    CompileParcellations(uiDCParameters,ListBox,AnalysisNameLabel,ParcellationLabel,VariableLabel,SelectNameButton,SelectParcelButton,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat,DataType) );
SelectTypeButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','Off','Text','Select Type','ButtonPushedFcn',@(SelectTypeButton,event) ...
    CompileAnalysisNames(uiDCParameters,ListBox,AnalysisTypeLabel,SelectTypeButton,SelectNameButton) );
CompileTypesButton=uibutton(UIFig,'Position',[20,10,100,40],'Visible','On','Text','Compile Types','ButtonPushedFcn',@(CompileTypesButton,event) ...
    CompileAnalysisTypes(uiDCParameters,ListBox,CompileTypesButton,SelectTypeButton,SubjectOrRun) );
waitfor(UIFig,'Visible','Off');

if uiDCParameters.ndim==1    
    [compiledData,DataInfo.compiled_fmriprep_table,DataInfo.selectind,DataInfo.DataLabels,DataInfo.ExperimentsDir,DataInfo.fmriprep_table,DataInfo.fmriprep_table_name] = ...
        Compile1D(uiDCParameters.ExperimentsDir,uiDCParameters.fmriprep_table,...
        'AnalysisType',uiDCParameters.AnalysisType,...
        'AnalysisName',uiDCParameters.AnalysisName,...
        'LoadVarName',uiDCParameters.LoadVarName,...
        'LoadVarFormat',uiDCParameters.LoadVarFormat,...
        'DataType',uiDCParameters.DataType,...
        'SubjectOrRun',uiDCParameters.SubjectOrRun,...
        'ParcelName',uiDCParameters.ParcelName,...
        'TableOrArray',uiDCParameters.TableOrArray);
elseif uiDCParameters.ndim==2
    [compiledData,DataInfo.compiled_fmriprep_table,DataInfo.selectind,DataInfo.DataLabels,DataInfo.ExperimentsDir,DataInfo.fmriprep_table,DataInfo.fmriprep_table_name] = ...
        Compile2D(uiDCParameters.ExperimentsDir,uiDCParameters.fmriprep_table,...
        'AnalysisType',uiDCParameters.AnalysisType,...
        'AnalysisName',uiDCParameters.AnalysisName,...
        'LoadVarName',uiDCParameters.LoadVarName,...
        'LoadVarFormat',uiDCParameters.LoadVarFormat,...
        'DataType',uiDCParameters.DataType,...
        'SubjectOrRun',uiDCParameters.SubjectOrRun,...
        'ParcelName',uiDCParameters.ParcelName,...
        'TableOrArray',uiDCParameters.TableOrArray,...
        'VertOrMat',uiDCParameters.VertOrMat); 
else
    [compiledData,DataInfo.compiled_fmriprep_table,DataInfo.selectind,DataInfo.DataLabels,DataInfo.ExperimentsDir,DataInfo.fmriprep_table,DataInfo.fmriprep_table_name] = ...
        CompileND(uiDCParameters.ExperimentsDir,uiDCParameters.fmriprep_table,...
        'AnalysisType',uiDCParameters.AnalysisType,...
        'AnalysisName',uiDCParameters.AnalysisName,...
        'LoadVarName',uiDCParameters.LoadVarName,...
        'LoadVarFormat',uiDCParameters.LoadVarFormat,...
        'DataType',uiDCParameters.DataType,...
        'SubjectOrRun',uiDCParameters.SubjectOrRun,...
        'ParcelName',uiDCParameters.ParcelName);     
end
clear global;
close all force;

end


function GetExpDirPush(uiDCParameters,ExpDir)
    ExperimentsDir=uigetdir('/','Select Experiments Directory');
    ExperimentsDir=[ExperimentsDir,'/'];
    ExpDir.Text=ExperimentsDir;
    uiDCParameters.ExperimentsDir=ExperimentsDir;
end

function GetTableNamePush(uiDCParameters,FmriPrepTableName)
    [~,fmriprep_table_names]=getFolderAndFileNames('fmriprep_table/');
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect([{'none- make table'};fmriprep_table_names],'Select fmriprep_table used:',SingleSelect,1);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
    FmriPrepTableName.Text=fmriprep_table_name;
    load(['fmriprep_table/',fmriprep_table_name],'fmriprep_table');
    uiDCParameters.fmriprep_table=fmriprep_table;
    uiDCParameters.fmriprep_table_name=fmriprep_table_name; 
end

function CompileAnalysisTypes(uiDCParameters,ListBox,CompileTypesButton,SelectTypeButton,SubjectOrRun)
    uiDCParameters.SubjectOrRun=SubjectOrRun.Value;
    ExperimentsDir=uiDCParameters.ExperimentsDir;
    fmriprep_table=uiDCParameters.fmriprep_table;
    TotalRuns=height(fmriprep_table);
    dataInd=find(fmriprep_table.run==1)';
    if strcmpi(uiDCParameters.SubjectOrRun,'Subject')
        useIndicies=dataInd;
    else
        useIndicies=[1:TotalRuns];
    end
    AllAnalysisTypes=[];
    for i = useIndicies
        dirName=[ExperimentsDir,fmriprep_table.matDir{i,1}];
        [folderNames]=getFolderAndFileNames(dirName);
        AllAnalysisTypes=unique([AllAnalysisTypes;folderNames]);
    end
    ListBox.Items=AllAnalysisTypes;
    CompileTypesButton.Visible='Off';
    SelectTypeButton.Visible='On';
    SubjectOrRun.Visible='Off';
end

function CompileAnalysisNames(uiDCParameters,ListBox,AnalysisTypeLabel,SelectTypeButton,SelectNameButton)
    uiDCParameters.AnalysisType=ListBox.Value;
    AnalysisTypeLabel.Text=[AnalysisTypeLabel.Text,uiDCParameters.AnalysisType];
    ExperimentsDir=uiDCParameters.ExperimentsDir;
    fmriprep_table=uiDCParameters.fmriprep_table;    
    TotalRuns=height(uiDCParameters.fmriprep_table);
    dataInd=find(fmriprep_table.run==1)';
    if strcmpi(uiDCParameters.SubjectOrRun,'Subject')
        useIndicies=dataInd;
    else
        useIndicies=[1:TotalRuns];
    end
    AllAnalysisNames=[];
    for i = useIndicies
        dirName=[ExperimentsDir,fmriprep_table.matDir{i,1},'/',uiDCParameters.AnalysisType,'/'];
        [folderNames]=getFolderAndFileNames(dirName);
        AllAnalysisNames=unique([AllAnalysisNames;folderNames]);
    end

    ListBox.Items=AllAnalysisNames;
    SelectTypeButton.Visible='Off';
    SelectNameButton.Visible='On';

end

function CompileParcellations(uiDCParameters,ListBox,AnalysisNameLabel,ParcellationLabel,VariableLabel,SelectNameButton,SelectParcelButton,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat,DataType)
    uiDCParameters.AnalysisName=ListBox.Value;
    AnalysisNameLabel.Text=[AnalysisNameLabel.Text,uiDCParameters.AnalysisName];
    [uiDCParameters.filePaths] = BIDsDirSearch(uiDCParameters.ExperimentsDir,uiDCParameters.fmriprep_table,'SubjectOrRun',uiDCParameters.SubjectOrRun,'AnalysisType',uiDCParameters.AnalysisType,'AnalysisName',uiDCParameters.AnalysisName);

    ListBox.Items=uiDCParameters.filePaths.Properties.VariableNames;
    SelectNameButton.Visible='Off';
    if length(ListBox.Items)>1
        SelectParcelButton.Visible='On';
        uiDCParameters.DataType='ByParcellation';
    else
        DataType.Visible='On';
        SelectVarButton.Visible='On';
        CompileVariables(uiDCParameters,ListBox,ParcellationLabel,VariableLabel,SelectParcelButton,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat);      
    end
end


function CompileVariables(uiDCParameters,ListBox,ParcellationLabel,VariableLabel,SelectParcelButton,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat)
    uiDCParameters.ParcelName=ListBox.Value;
    ParcellationLabel.Text=[ParcellationLabel.Text,uiDCParameters.ParcelName];
    uiDCParameters.filePaths=uiDCParameters.filePaths.(uiDCParameters.ParcelName);  
    for i = 1:size(uiDCParameters.filePaths,1)
        variableInfo = who('-file', uiDCParameters.filePaths{i,1});
        LoadVariableName=unique([LoadVariableName;variableInfo]);
    end
    ListBox.Items=LoadVariableName;
    SelectParcelButton.Visible='Off';
    if length(ListBox.Items)>1
        SelectVarButton.Visible='On';
    else
        CompileVarInfo(uiDCParameters,ListBox,VariableLabel,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat)
    end
end

function CompileVarInfo(uiDCParameters,ListBox,VariableLabel,SelectVarButton,CompileDataButton,TableOrArray,VertOrMat)
    SelectVarButton.Visible='Off';    
    CompileDataButton.Visible='On';
    uiDCParameters.VarName=ListBox.Value;
    VariableLabel.Text=[VariableLabel.Text,uiDCParameters.VarName];
   
    for i = 1:size(uiDCParameters.filePaths,1)
        try
            tempLoad=load(uiDCParameters.filePaths{i,1},uiDCParameters.VarName);
            tempVar=tempLoad.(uiDCParameters.VarName);
            if istable(tempVar)
                uiDCParameters.LoadVarFormat='Table';
                TableOrArray.Visible='On';
                if any(ismember(size(tempVar),1))
                    uiDCParameters.ndim=1;
                else
                    uiDCParameters.ndim=2;
                    VertOrMat.Visible='On';
                end
            elseif iscell(tempVar)
                uiDCParameters.LoadVarFormat='Cell';
                if ndims(tempVar)>2
                    uiDCParameters.ndim=ndims(tempVar);
                elseif any(ismember(size(tempVar),1))
                    uiDCParameters.ndim=1;
                else  
                    uiDCParameters.ndim=2;
                    VertOrMat.Visible='On';
                end
            else
                uiDCParameters.LoadVarFormat='Array';
                if ndims(tempVar)>2
                    uiDCParameters.ndim=ndims(tempVar);
                elseif any(ismember(size(tempVar),1))
                    uiDCParameters.ndim=1;
                else  
                    uiDCParameters.ndim=2;
                    VertOrMat.Visible='On';
                end
            end
            break
        catch
            continue
        end
    end 
end

function CompileDataPush(uiDCParameters,UIFig,VertOrMat,TableOrArray,DataType)
    uiDCParameters.DataType=DataType.Value;
    uiDCParameters.VertOrMat=VertOrMat.Value;
    uiDCParameters.TableOrArray=TableOrArray.Value;
    UIFig.Visible='Off';
end