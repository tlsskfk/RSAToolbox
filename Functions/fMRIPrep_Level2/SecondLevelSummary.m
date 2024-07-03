function [Results,ResultsByType,ResultsFullFeature] = SecondLevelSummary(AllTables,ParcelNames)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numParcels=size(AllTables,1);
numParcelsOrig=numParcels;
NodeVarNamesAll=[];
EdgeVarNamesAll=[];

TType=1:3;
for t = TType
    for i = 1:numParcels
        if ~isempty(AllTables{i,t})
            NodeVarNamesAll=unique([NodeVarNamesAll;AllTables{i,t}.Properties.VariableNames']);
        end
        try
        if ~isempty(AllTables{i,t+3})
            EdgeVarNamesAll=unique([EdgeVarNamesAll;AllTables{i,t+3}.Properties.VariableNames']);
        end
        end
    end
end
UseVarNames=[];
for t = TType    
    for i = 1:numParcels
       if ~isempty(AllTables{i,t})
           [~,tempVarNames] = StringMatch(AllTables{i,t}.Properties.VariableNames','_');
%            tempVarNames=AllTables{i,t}.Properties.VariableNames';
            if isempty(UseVarNames)
                UseVarNames=tempVarNames;
            else
                UseVarNames=tempVarNames(ismember(tempVarNames,UseVarNames));
            end
       end
       try
        if ~isempty(AllTables{i,t+3})
            [~,tempVarNames] = StringMatch(AllTables{i,t+3}.Properties.VariableNames','_');
            if isempty(UseVarNames)
                UseVarNames=tempVarNames;
            else
                UseVarNames=tempVarNames(ismember(tempVarNames,UseVarNames));
            end
        end    
       end
    end
end
[Prefix] = StringMatch([NodeVarNamesAll(:);EdgeVarNamesAll(:)],'_');
NodePrefix=uiNameSelect(Prefix,'Select node prefixes.');
EdgePrefix=uiNameSelect(Prefix,'Select edge prefixes.');
VarNames=uiNameSelect(UseVarNames,'Select variable names.');
NodeVarNames=[];
EdgeVarNames=[];
for i = 1:length(NodePrefix)
    tempCat=[repmat(NodePrefix(i,1),[length(VarNames),1]),VarNames];
    NodeVarNames=[NodeVarNames;join(tempCat,'')];
end
for i = 1:length(EdgePrefix)
    tempCat=[repmat(EdgePrefix(i,1),[length(VarNames),1]),VarNames];
    EdgeVarNames=[EdgeVarNames;join(tempCat,'')];
end



Results=cell(1);
ResultsByType=cell(1);
ResultsFullFeature=cell(1);
for t = TType
    numParcels=numParcelsOrig;
    count=1;
    AllNodeTable=[];
    AllEdgeTable=[];
    for i = 1:numParcels
        NodeTable=AllTables{i,t};
        NodeTable.Properties.RowNames=join([NodeTable.Properties.RowNames,repmat(ParcelNames(i,1),[length(NodeTable.Properties.RowNames),1])],'_');
        AllNodeTable=[AllNodeTable;NodeTable];
        
        EdgeTable=AllTables{i,t+3};
        EdgeTable.Properties.RowNames=join([EdgeTable.Properties.RowNames,repmat(ParcelNames(i,1),[length(EdgeTable.Properties.RowNames),1])],'_');
        AllEdgeTable=[AllEdgeTable;EdgeTable];
    end
    if count==1
        numParcels=numParcels+1;
    end
    count=count+1;
    AllTables{numParcels,t}=AllNodeTable;
    AllTables{numParcels,t+3}=AllEdgeTable;   
    NodeTableByType=cell(numParcels,length(NodePrefix));
    EdgeTableByType=cell(numParcels,length(EdgePrefix));
    for i = 1:numParcels
        AllFeatureTable=[];
        AllFeatureTable_Node=[];
        AllFeatureTable_Edge=[];
        NodeTable=AllTables{i,t};
        EdgeTable=AllTables{i,t+3};
        UseNodeVars=NodeVarNames(ismember(NodeVarNames,NodeTable.Properties.VariableNames));
        UseEdgeVars=EdgeVarNames(ismember(EdgeVarNames,EdgeTable.Properties.VariableNames));
        NodeTable=NodeTable(:,UseNodeVars);
        EdgeTable=EdgeTable(:,UseEdgeVars);
        try
            [Results{i,t}.Table,Results{i,t}.TableT,Results{i,t}.ClusterResults]=ComputeTableMarginalSums(NodeTable);
        catch
            disp('error')
        end
        try
            [Results{i,t+3}.Table,Results{i,t+3}.TableT,Results{i,t+3}.ClusterResults]=ComputeTableMarginalSums(EdgeTable);
        catch
            disp('error')
        end
        for j = 1:length(NodePrefix)
            tempNodeVars=UseNodeVars(contains(UseNodeVars,NodePrefix{j,1}));
            tempNodeTable=NodeTable(:,tempNodeVars);
            tempNodeTable.Properties.VariableNames=strrepCell(tempNodeVars,NodePrefix{j,1},'');
            tempNodeTable.Properties.RowNames=join([repmat(NodePrefix(j,1),[length(tempNodeTable.Properties.RowNames),1]),tempNodeTable.Properties.RowNames],'');
            AllFeatureTable=[AllFeatureTable;tempNodeTable];
            AllFeatureTable_Node=[AllFeatureTable_Node;tempNodeTable];
            try
                [ResultsByType{i,t}.(NodePrefix{j,1}).Table,ResultsByType{i,t}.(NodePrefix{j,1}).TableT,ResultsByType{i,t}.(NodePrefix{j,1}).ClusterResults]=ComputeTableMarginalSums(tempNodeTable);
            catch
                disp('error')
            end
        end
        for j = 1:length(EdgePrefix)
            tempEdgeVars=UseEdgeVars(contains(UseEdgeVars,EdgePrefix{j,1}));
            tempEdgeTable=EdgeTable(:,tempEdgeVars);
            tempEdgeTable.Properties.VariableNames=strrepCell(tempEdgeVars,EdgePrefix{j,1},'');
            tempEdgeTable.Properties.RowNames=join([repmat(EdgePrefix(j,1),[length(tempEdgeTable.Properties.RowNames),1]),tempEdgeTable.Properties.RowNames],'');
            AllFeatureTable=[AllFeatureTable;tempEdgeTable];
            AllFeatureTable_Edge=[AllFeatureTable_Edge;tempEdgeTable];            
            try
                [ResultsByType{i,t+3}.(EdgePrefix{j,1}).Table,ResultsByType{i,t+3}.(EdgePrefix{j,1}).TableT,ResultsByType{i,t+3}.(EdgePrefix{j,1}).ClusterResults]=ComputeTableMarginalSums(tempEdgeTable);
            catch
                disp('error')
            end
        end  
        try
            [ResultsFullFeature{i,t}.Node.Table,ResultsFullFeature{i,t}.Node.TableT,ResultsFullFeature{i,t}.Node.ClusterResults]=ComputeTableMarginalSums(AllFeatureTable_Node);
        catch
            disp('error')
        end
        try
            [ResultsFullFeature{i,t}.Edge.Table,ResultsFullFeature{i,t}.Edge.TableT,ResultsFullFeature{i,t}.Edge.ClusterResults]=ComputeTableMarginalSums(AllFeatureTable_Edge);
        catch
            disp('error')
        end
        try
            [ResultsFullFeature{i,t}.All.Table,ResultsFullFeature{i,t}.All.TableT,ResultsFullFeature{i,t}.All.ClusterResults]=ComputeTableMarginalSums(AllFeatureTable);
        catch
            disp('error')
        end
    end 
end

end

