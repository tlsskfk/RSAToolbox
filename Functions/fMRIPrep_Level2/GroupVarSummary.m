function [Results,ParcelNames] = GroupVarSummary(ExperimentsDir,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[editTableNames] = VariableSetter('editTableNames',1,varargin);
[NoSave] = VariableSetter('NoSave',0,varargin);
SaveDir=[ExperimentsDir,'Group_FeatureMaps/'];

[AllPaths,~,AllFolderNames,~,AllInfo] = GroupDirSearch(ExperimentsDir);
UseVarNames=AllInfo.LoadVarNames_Parcel;
NodeVarNames=uiNameSelect([{'None'};UseVarNames],'Select node variable names:',0);
if strcmpi(NodeVarNames,'None')
    NodeVarNames=[];
else
    UseVarNames=UseVarNames(~ismember(UseVarNames,NodeVarNames),:);
end
EdgeVarNames=uiNameSelect([{'None'};UseVarNames],'Select edge variable names:',0);        
ParcelNames=AllInfo.ParcellationNames;

Paths_Nodes=[];
FolderNames_Nodes=[];
Paths_Edges=[];
FolderNames_Edges=[];
for i = 1:height(AllFolderNames)
    if ~isempty(NodeVarNames)
        if any(ismember(AllFolderNames.VarNames{i,1},NodeVarNames))
            tempNodeVars=AllFolderNames.VarNames{i,1}(ismember(AllFolderNames.VarNames{i,1},NodeVarNames),:);
            if ~iscell(tempNodeVars)
                tempNodeVars={tempNodeVars};
            end
            for j = 1:size(tempNodeVars,1)
                tempPaths=AllPaths(i,:);
                tempPaths{1,2}=tempNodeVars{j,1};
                Paths_Nodes=[Paths_Nodes;tempPaths];
                tempFolderNames=AllFolderNames(i,:);
                tempFolderNames.VarNames{1,1}=tempNodeVars{j,1};
                FolderNames_Nodes=[FolderNames_Nodes;tempFolderNames];  
            end
        end
    end
    
    if ~isempty(EdgeVarNames)
        if any(ismember(AllFolderNames.VarNames{i,1},EdgeVarNames))
            tempEdgeVars=AllFolderNames.VarNames{i,1}(ismember(AllFolderNames.VarNames{i,1},EdgeVarNames),:);
            if ~iscell(tempEdgeVars)
                tempEdgeVars={tempEdgeVars};
            end
            for j = 1:size(tempEdgeVars,1)
                tempPaths=AllPaths(i,:);
                tempPaths{1,2}=tempEdgeVars{j,1};
                Paths_Edges=[Paths_Edges;tempPaths];
                tempFolderNames=AllFolderNames(i,:);
                tempFolderNames.VarNames{1,1}=tempEdgeVars{j,1};
                FolderNames_Edges=[FolderNames_Edges;tempFolderNames];  
            end
        end
    end
end

UseIDVarNames_Nodes=[];
UseIDVarNames_Edges=[];
FigureFolderName_Nodes=[];
FigureFileName_Nodes=[];
FigureFolderName_Edges=[];
FigureFileName_Edges=[];
UniqueNames_Nodes=[];
for i = [1:3,5,6]
    if ~isempty(Paths_Nodes)
        UniqueNames_Nodes=table2cell(unique(FolderNames_Nodes(:,i)));
        UniqueNames_Nodes(strcmpi(UniqueNames_Nodes,'AnalysisParameters.mat'),:)=[];
        if length(UniqueNames_Nodes)>1
            tempTableNodes=FolderNames_Nodes.(FolderNames_Nodes.Properties.VariableNames{i});
            tempTableNodes=strrepCell(tempTableNodes,'.mat','');
            if editTableNames==1
                [VarOrder_Nodes,NewVarNames_Nodes] = uiTableVarNames(UniqueNames_Nodes);
                VarOrder_Nodes=VarOrder_Nodes(:);
                NewVarNames_Nodes=NewVarNames_Nodes(:);
                for j = 1: size(NewVarNames_Nodes,1)
                    tempTableNodes=strrepCell(tempTableNodes,VarOrder_Nodes{j,1},NewVarNames_Nodes{j,1});
                end
            end
            UseIDVarNames_Nodes=[UseIDVarNames_Nodes,tempTableNodes];
            FigureFolderName_Nodes=[FigureFolderName_Nodes,FolderNames_Nodes.Properties.VariableNames{i},'_'];
        else
            FigureFileName_Nodes=[FigureFileName_Nodes,FolderNames_Nodes.(FolderNames_Nodes.Properties.VariableNames{i}){1,1},'_'];
        end
    end
    if ~isempty(Paths_Edges)
        UniqueNames_Edges=table2cell(unique(FolderNames_Edges(:,i)));
        UniqueNames_Edges(strcmpi(UniqueNames_Edges,'AnalysisParameters.mat'),:)=[];
        if length(UniqueNames_Edges)>1
            tempTableEdges=FolderNames_Edges.(FolderNames_Edges.Properties.VariableNames{i});
            tempTableEdges=strrepCell(tempTableEdges,'.mat','');
            if editTableNames==1
                if ~isempty(UniqueNames_Nodes) && size(UniqueNames_Nodes,1) == size(UniqueNames_Edges,1) && all(ismember(UniqueNames_Nodes,UniqueNames_Edges))
                    NewVarNames_Edges=NewVarNames_Nodes;
                    VarOrder_Edges=VarOrder_Nodes;
                else
                    [VarOrder_Edges,NewVarNames_Edges] = uiTableVarNames(UniqueNames_Edges);
                    VarOrder_Edges=VarOrder_Edges(:);
                    NewVarNames_Edges=NewVarNames_Edges(:);
                end
                for j = 1: size(NewVarNames_Edges,1)
                    tempTableEdges=strrepCell(tempTableEdges,VarOrder_Edges{j,1},NewVarNames_Edges{j,1});
                end
            end
            UseIDVarNames_Edges=[UseIDVarNames_Edges,tempTableEdges];
            FigureFolderName_Edges=[FigureFolderName_Edges,FolderNames_Edges.Properties.VariableNames{i},'_'];
        else
            FigureFileName_Edges=[FigureFileName_Edges,FolderNames_Edges.(FolderNames_Edges.Properties.VariableNames{i}){1,1},'_'];
        end   
    end
end

if ~isempty(Paths_Nodes)
    if isempty(FigureFolderName_Nodes)
        FigureFolderName_Nodes='SingleGroup_Nodes';
    end
    if isempty(FigureFileName_Nodes) && NoSave==0
        FigureFileName_Nodes=uiEnterName('','Enter node filename');
    end    
    if isempty(UseIDVarNames_Nodes)
        UseIDVarNames_Nodes{1,1}=uiEnterName(FigureFileName_Nodes,'Enter node var name');
        NewVarNames_Nodes=UseIDVarNames_Nodes;
    end
    if NoSave==0
        FigureFolderName_Nodes(:,end)='';
        FigureFileName_Nodes(:,end)='';
        FigureFileName_Nodes=strrep(FigureFileName_Nodes,'.mat','');
    end
    if size(UseIDVarNames_Nodes,2)>1
        UseIDVarNames_Nodes=join(UseIDVarNames_Nodes,2);
        UseIDVarNames_Nodes=strrepCell(UseIDVarNames_Nodes,' ','_');
    end
    pVars_Nodes=Paths_Nodes(:,2);
    statVars_Nodes=Paths_Nodes(:,2);
    for i = 1:size(NodeVarNames,1)
        tempPaths_Nodes=Paths_Nodes(strcmp(Paths_Nodes(:,2),NodeVarNames{i,1}),:);
        tempTable=load(tempPaths_Nodes{1,1},Paths_Nodes{1,2});
        tempTable=tempTable.(Paths_Nodes{1,2});
        tempVarNames=tempTable.Properties.VariableNames;
        pVars_Name=uiNameSelect(tempVarNames(:),['Select node p value for ',Paths_Nodes{1,2}],1);
        pVars_Nodes=strrepCell(pVars_Nodes,Paths_Nodes{1,2},pVars_Name);
        statVars_Name=uiNameSelect(tempVarNames(:),['Select node stat value for ',Paths_Nodes{1,2}],1);
        statVars_Nodes=strrepCell(statVars_Nodes,Paths_Nodes{1,2},statVars_Name);
    end
    Paths_Nodes=[Paths_Nodes,pVars_Nodes,statVars_Nodes];
    if NoSave==0
        AnalysisName_Nodes=uiEnterName(FigureFileName_Nodes,'Enter nodes analysis name:');
    end
end

if ~isempty(Paths_Edges)
    if isempty(FigureFolderName_Edges)
        FigureFolderName_Edges='SingleGroup_Edges';
    end
    if isempty(FigureFileName_Edges) && NoSave==0
        FigureFileName_Edges=uiEnterName('','Enter edge filename');
    end    
    if isempty(UseIDVarNames_Edges)
        UseIDVarNames_Edges{1,1}=uiEnterName(FigureFileName_Edges,'Enter edge var name');
        NewVarNames_Edges=UseIDVarNames_Edges;
    end
    if NoSave==0
        FigureFolderName_Edges(:,end)='';
        FigureFileName_Edges(:,end)='';
        FigureFileName_Edges=strrep(FigureFileName_Edges,'.mat','');
    end    
    if size(UseIDVarNames_Edges,2)>1
        UseIDVarNames_Edges=join(UseIDVarNames_Edges,2);
        UseIDVarNames_Edges=strrepCell(UseIDVarNames_Edges,' ','_');
    end
    pVars_Edges=Paths_Edges(:,2);
    statVars_Edges=Paths_Edges(:,2);
    for i = 1:size(EdgeVarNames,1)
        tempPaths_Edges=Paths_Edges(strcmp(Paths_Edges(:,2),EdgeVarNames{i,1}),:);
        tempTable=load(tempPaths_Edges{1,1},Paths_Edges{1,2});
        tempTable=tempTable.(Paths_Edges{1,2});
        tempVarNames=tempTable.Properties.VariableNames;
        pVars_Name=uiNameSelect(tempVarNames(:),['Select edge p value for ',Paths_Edges{1,2}],1);
        pVars_Edges=strrepCell(pVars_Edges,Paths_Edges{1,2},pVars_Name);
        statVars_Name=uiNameSelect(tempVarNames(:),['Select edge stat value for ',Paths_Edges{1,2}],1);
        statVars_Edges=strrepCell(statVars_Edges,Paths_Edges{1,2},statVars_Name);
    end
    Paths_Edges=[Paths_Edges,pVars_Edges,statVars_Edges];
    if NoSave==0
        AnalysisName_Edges=uiEnterName(FigureFileName_Edges,'Enter edges analysis name:');
    end
end

GroupDir=strrep(ExperimentsDir,'/Experiments','');
Results=cell(length(ParcelNames),2);
for parcelNum=1:length(ParcelNames)
    ParcelName=ParcelNames{parcelNum,1};
    try
        load(['Parcellations/',ParcelName],'UseLabels');
    catch
        try
            [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
            load(['Parcellations/IndvParcels/',fmriprep_table_name,'/',ParcelName],'UseLabels');
        catch
            load(['Parcellations/IndvCoords/',fmriprep_table_name,'/',ParcelName],'UseLabels');
        end
    end
    
    if ~isempty(Paths_Nodes)
        
        SaveDir_Nodes=[GroupDir,'Group_FeatureMaps_Nodes/',FigureFolderName_Nodes,'/',ParcelName,'/'];
        if NoSave==0
            if ~exist(SaveDir_Nodes,'file')
                mkdir(SaveDir_Nodes);
            end        
        end      
        UsePaths_Nodes=Paths_Nodes(strcmp(FolderNames_Nodes.ParcellationNames,ParcelName),:); 
        NodeTable=array2table(zeros(length(UseLabels),size(UsePaths_Nodes,1)),'VariableNames',NewVarNames_Nodes,'RowNames',UseLabels);   
        NodeTableBonf=array2table(zeros(length(UseLabels),size(UsePaths_Nodes,1)),'VariableNames',NewVarNames_Nodes,'RowNames',UseLabels); 
        NodeTableRaw=array2table(zeros(length(UseLabels),size(UsePaths_Nodes,1)),'VariableNames',NewVarNames_Nodes,'RowNames',UseLabels);          
        for i = 1:size(UsePaths_Nodes,1)
            tempTable=load(UsePaths_Nodes{i,1},UsePaths_Nodes{i,2});
            [ SigMat ] = SigSeg2( tempTable.(UsePaths_Nodes{i,2}).(UsePaths_Nodes{i,3}),[0.1,0.05,0.01,0.005,0.001],'pos');
            NodeTable.(NewVarNames_Nodes{i,1})=SigMat.*sign(tempTable.(UsePaths_Nodes{i,2}).(UsePaths_Nodes{i,4}));
            [ SigMat ] = SigSeg2( tempTable.(UsePaths_Nodes{i,2}).(UsePaths_Nodes{i,3}),[0.05/length(UseLabels)],'pos');
            NodeTableBonf.(NewVarNames_Nodes{i,1})=SigMat.*sign(tempTable.(UsePaths_Nodes{i,2}).(UsePaths_Nodes{i,4})); 
            NodeTableRaw.(NewVarNames_Nodes{i,1})=tempTable.(UsePaths_Nodes{i,2}).(UsePaths_Nodes{i,4});            
        end 
        if NoSave==0
            save([SaveDir_Nodes,AnalysisName_Nodes],'NodeTable','NodeTableBonf','NodeTableRaw')
        end    
        Results{parcelNum,1}=NodeTable;
        Results{parcelNum,2}=NodeTableBonf;
        Results{parcelNum,3}=NodeTableRaw;
    end
    if ~isempty(Paths_Edges)
        SaveDir_Edges=[GroupDir,'Group_FeatureMaps_Edges/',FigureFolderName_Edges,'/',ParcelName,'/'];
        if NoSave==0
            if ~exist(SaveDir_Edges,'file')
                mkdir(SaveDir_Edges);
            end         
        end
        [ labelPairs1 ] = labels2uppertriuvectorlabels( UseLabels );     
        UsePaths_Edges=Paths_Edges(strcmp(FolderNames_Edges.ParcellationNames,ParcelName),:); 
        EdgeTable=array2table(zeros(length(labelPairs1),size(UsePaths_Edges,1)),'VariableNames',NewVarNames_Edges,'RowNames',labelPairs1);    
        EdgeTableBonf=array2table(zeros(length(labelPairs1),size(UsePaths_Edges,1)),'VariableNames',NewVarNames_Edges,'RowNames',labelPairs1);  
        EdgeTableRaw=array2table(zeros(length(labelPairs1),size(UsePaths_Edges,1)),'VariableNames',NewVarNames_Edges,'RowNames',labelPairs1);             
        for i = 1:size(UsePaths_Edges,1)
            tempTable=load(UsePaths_Edges{i,1},UsePaths_Edges{i,2});
            [ SigMat ] = SigSeg2( tempTable.(UsePaths_Edges{i,2}).(UsePaths_Edges{i,3}),[0.1,0.05,0.01,0.005,0.001],'pos');
            EdgeTable.(NewVarNames_Edges{i,1})=SigMat.*sign(tempTable.(UsePaths_Edges{i,2}).(UsePaths_Edges{i,4}));
            [ SigMat ] = SigSeg2( tempTable.(UsePaths_Edges{i,2}).(UsePaths_Edges{i,3}),[0.05/length(labelPairs1)],'pos');
            EdgeTableBonf.(NewVarNames_Edges{i,1})=SigMat.*sign(tempTable.(UsePaths_Edges{i,2}).(UsePaths_Edges{i,4}));   
            EdgeTableRaw.(NewVarNames_Edges{i,1})=tempTable.(UsePaths_Edges{i,2}).(UsePaths_Edges{i,4});
        end    
        if NoSave==0
            save([SaveDir_Edges,AnalysisName_Edges],'EdgeTable','EdgeTableBonf','EdgeTableRaw')
        end    
        Results{parcelNum,4}=EdgeTable;
        Results{parcelNum,5}=EdgeTableBonf;
        Results{parcelNum,6}=EdgeTableRaw;
    end    
end
end
