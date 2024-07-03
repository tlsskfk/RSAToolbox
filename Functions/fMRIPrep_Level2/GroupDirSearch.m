function [AllPaths_Parcel,AllPaths_NoParcel,FolderNames_Parcel,FolderNames_NoParcel,AllInfo] = GroupDirSearch(ExperimentsDir,varargin)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[TableNames] = VariableSetter('TableNames',[],varargin);
[AnalysisTypes] = VariableSetter('AnalysisTypes',[],varargin);
[AnalysisNames] = VariableSetter('AnalysisNames',[],varargin);
[ParcellationNames] = VariableSetter('ParcellationNames',[],varargin);
[LoadVarNames_NoParcel] = VariableSetter('LoadVarNames_NoParcel',[],varargin);
[LoadVarNames_Parcel] = VariableSetter('LoadVarNames_Parcel',[],varargin);
[FileNames] = VariableSetter('FileNames',[],varargin);
GroupType = VariableSetter('GroupType',[''],varargin);
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'/Experiments/','');
AllInfo=struct;
if isempty(GroupType)
    [GroupType] = uiNameSelect({'GroupAnalysis','GroupTables','IndvDiff_GroupAnalysis','GroupDiff_GroupAnalysis'},'Select type of group folder to search',1);
end
if iscell(GroupType)
    GroupType=GroupType{1,1};
end
GroupDir=[GroupDir,'/',GroupType,'/'];
GroupDir=strrep(GroupDir,'//','/');



if isempty(TableNames)
    TableNames=getFolderAndFileNames(GroupDir);
    SingleSelect=0;
    [TableNames] = uiNameSelect(TableNames,'Select table names for file list.',SingleSelect);
else
    if ~iscell(TableNames)
        TableNames={TableNames};
    end
    TableNames=TableNames(:);
end
        
if isempty(AnalysisTypes)
    SingleSelect=0;
    for i = 1:length(TableNames)
        dirName=[GroupDir,TableNames{i,1},'/'];
        [folderNames]=getFolderAndFileNames(dirName);
        AnalysisTypes=unique([AnalysisTypes;folderNames]);
    end
    [AnalysisTypes] = uiNameSelect(AnalysisTypes,'Select analysis types to include:' ,SingleSelect);
else
    if ~iscell(AnalysisTypes)
        AnalysisTypes={AnalysisTypes};
    end
    AnalysisTypes=AnalysisTypes(:);
end
[AllDirNames] = FolderNames2Dirs({GroupDir,TableNames,AnalysisTypes});
for i = 1:length(AllDirNames)
    if ~exist(AllDirNames{i,1},'dir')
        AllDirNames(i,:)=[];
    end
end

if isempty(AnalysisNames)
    SingleSelect=0;
    for i = 1:length(AllDirNames)
        dirName=[AllDirNames{i,1}];
        [folderNames]=getFolderAndFileNames(dirName);
        AnalysisNames=unique([AnalysisNames;folderNames]);
    end
    [AnalysisNames] = uiNameSelect(AnalysisNames,'Select analysis names to include:' ,SingleSelect);
else
    if ~iscell(AnalysisNames)
        AnalysisNames={AnalysisNames};
    end
    AnalysisNames=AnalysisNames(:);
end
[AllDirNames] = FolderNames2Dirs({AllDirNames,AnalysisNames});
temRemInd=zeros(length(AllDirNames),1);
for i = 1:length(AllDirNames)
    if ~exist(AllDirNames{i,1},'dir')
        temRemInd(i,1)=1;
    end
end
AllDirNames(temRemInd==1,:)=[];
tempFolderNames=[];
AllFileNames=[];
for i = 1:length(AllDirNames)
    dirName=[AllDirNames{i,1}];
    [folderNames,tempfileNames]=getFolderAndFileNames(dirName);
    tempFolderNames=unique([tempFolderNames;folderNames]);
    AllFileNames=unique([AllFileNames;tempfileNames]);
end

AllPaths_NoParcel=[];
if ~isempty(AllFileNames)
    SingleSelect=0;
    if isempty(FileNames)
        [AllFileNames] = uiNameSelect(AllFileNames,'Select files to include:',SingleSelect);
    else
        AllFileNames=FileNames;
    end
    [AllPaths_NoParcel] = FolderNames2Dirs({AllDirNames},AllFileNames);
end
for i = 1:length(AllPaths_NoParcel)
    if ~exist(AllPaths_NoParcel{i,1},'file')
        AllPaths_NoParcel(i,:)=[];
    end
end    
    
if ~isempty(tempFolderNames)    
    if isempty(ParcellationNames)
        SingleSelect=0;
        for i = 1:length(AllDirNames)
            dirName=[AllDirNames{i,1}];
            [folderNames]=getFolderAndFileNames(dirName);
            ParcellationNames=unique([ParcellationNames;folderNames]);
        end
        [ParcellationNames] = uiNameSelect(ParcellationNames,'Select parcellations to include:' ,SingleSelect);
    else
        if ~iscell(ParcellationNames)
            ParcellationNames={ParcellationNames};
        end
        ParcellationNames=ParcellationNames(:);
    end
    [AllDirNames] = FolderNames2Dirs({AllDirNames,ParcellationNames});
end
removeInd=zeros(length(AllDirNames),1);
for i = 1:length(AllDirNames)
    if ~exist(AllDirNames{i,1},'dir')
        removeInd(i,1)=1;
    end
end
AllDirNames(removeInd==1,:)=[];
AllFileNames=[];

for i = 1:length(AllDirNames)
    dirName=[AllDirNames{i,1}];
    [~,tempfileNames]=getFolderAndFileNames(dirName);
    AllFileNames=unique([AllFileNames;tempfileNames]);
end
if ~isempty(AllFileNames)
    SingleSelect=0;
    if isempty(FileNames)
        [AllFileNames] = uiNameSelect(AllFileNames,'Select files to include:',SingleSelect);
    else
        AllFileNames=FileNames;
    end
    [AllPaths_Parcel] = FolderNames2Dirs({AllDirNames},AllFileNames);
end
RemovePath=zeros(length(AllPaths_Parcel),1);
for i = 1:length(AllPaths_Parcel)
    if ~exist(AllPaths_Parcel{i,1},'file')
        RemovePath(i,1)=1;
    end
end   
AllPaths_Parcel(RemovePath==1,:)=[];
if ~isempty(AllPaths_Parcel)
	tempPaths = strrepCell(AllPaths_Parcel,GroupDir,'');
    FolderNames_Parcel = Dirs2FolderNames(tempPaths);
    FolderNames_Parcel = cell2table(FolderNames_Parcel,'VariableNames',{'TableNames','AnalysisTypes','AnalysisNames','ParcellationNames','FileNames'});
else
    FolderNames_Parcel=[];
end

if ~isempty(AllPaths_Parcel)
    VarNames_Parcel=cell(length(AllPaths_Parcel),1);    
    if isempty(LoadVarNames_Parcel)
        for i = 1:length(AllPaths_Parcel)        
            variableNames = who('-file', AllPaths_Parcel{i,1});
            VarNames_Parcel{i,1}=variableNames(:);
            LoadVarNames_Parcel=unique([LoadVarNames_Parcel;variableNames(:)]);
        end
        [LoadVarNames_Parcel] = uiNameSelect(LoadVarNames_Parcel,'Select parcellated variables to include:',SingleSelect);
    end
    for i = 1:length(AllPaths_Parcel)  
        VarNames_Parcel{i,1}=VarNames_Parcel{i,1}(ismember(VarNames_Parcel{i,1},LoadVarNames_Parcel),:);
    end
    FolderNames_Parcel=[FolderNames_Parcel,cell2table(VarNames_Parcel,'VariableNames',{'VarNames'})];
    AllPaths_Parcel=[AllPaths_Parcel,VarNames_Parcel];
end

if ~isempty(AllPaths_NoParcel)
	tempPaths = strrepCell(AllPaths_NoParcel,GroupDir,'');
    FolderNames_NoParcel = Dirs2FolderNames(tempPaths);
    FolderNames_NoParcel = cell2table(FolderNames_NoParcel,'VariableNames',{'TableNames','AnalysisTypes','AnalysisNames','FileNames'});
else
    FolderNames_NoParcel=[];
end

if ~isempty(AllPaths_NoParcel)
    VarNames_NoParcel=cell(length(AllPaths_NoParcel),1);
    if isempty(LoadVarNames_NoParcel)
        for i = 1:length(AllPaths_NoParcel)        
            variableNames = who('-file', AllPaths_NoParcel{i,1});
            VarNames_NoParcel{i,1}=variableNames(:);
            LoadVarNames_NoParcel=unique([LoadVarNames_NoParcel;variableNames(:)]);
        end
        [LoadVarNames_NoParcel] = uiNameSelect(LoadVarNames_NoParcel,'Select non parcellated variables to include:',SingleSelect);
    end
    for i = 1:length(AllPaths_Parcel)  
        VarNames_NoParcel{i,1}=VarNames_NoParcel{i,1}(ismember(VarNames_NoParcel{i,1},LoadVarNames_NoParcel),:);
    end
    FolderNames_NoParcel=[FolderNames_NoParcel,cell2table(VarNames_NoParcel,'VariableNames',{'VarNames'})];
    AllPaths_NoParcel=[AllPaths_Parcel,LoadVarNames_NoParcel];
end
FileNames=AllFileNames;

AllInfo.GroupDir=GroupDir;
AllInfo.AnalysisNames=AnalysisNames;
AllInfo.AnalysisTypes=AnalysisTypes;
AllInfo.TableNames=TableNames;
AllInfo.FileNames=FileNames;
AllInfo.ParcellationNames=ParcellationNames;
AllInfo.GroupType=GroupType;
AllInfo.LoadVarNames_Parcel=LoadVarNames_Parcel;
AllInfo.LoadVarNames_NoParcel=LoadVarNames_NoParcel;
end

