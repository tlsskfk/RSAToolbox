function [ Clusters ] = ClusterFinder(vals,maskcoords,posneg,threshtype,thresholds,clustersizethreshold,clusterdef) %'side';'line';'point 
%Written by David Rothlein
%NOTE: ZSCORE THRTESHOLDS  .1=1.645, .05=1.96, .01=2.576, .005=2.807, .001=3.291  
if strcmp(thresholds,'basic')==1
    thresholds=[1.96;2.576;2.807;3.291];
end;

if strcmp(threshtype,'count')==1

    Clusters=struct();
    Clusters.ThresholdedClusters=cell(1,1);
    numofclusts=max(vals);
    for i = 1:numofclusts
        tempClusters{i,1}=maskcoords(vals==i,:);
        tempClusters{i,2}=ones(size(tempClusters{i,1},1),1)*i;
        tempClusters{i,3}=size(tempClusters{i,1},1);
    end;
    Clusters.ThresholdedClusters{1,1} = tempClusters;
else
    
if size(thresholds,2)>1
    thresholds=thresholds';
end;

Clusters=struct();
Clusters.ThresholdedClusters=cell(1,length(thresholds));
Clusters.AllThresholdedCoords=cell(1,length(thresholds));
Clusters.MaxClusterSizes=zeros(1,length(thresholds));
Clusters.AllClusterSizes=cell(1,length(thresholds));

for i = 1:length(thresholds)

    if strcmp(posneg,'pos')==1
        [TempClusters,UncorrectedThresholdCoords,clustersizes,allcoords,allvals]=ClusterIdentifierSimple(vals,maskcoords,threshtype,thresholds(i,1),clustersizethreshold,clusterdef);
    elseif strcmp(posneg,'neg') == 1
        [TempClusters,UncorrectedThresholdCoords,clustersizes,allcoords,allvals]=ClusterIdentifierSimple(vals*-1,maskcoords,threshtype,thresholds(i,1),clustersizethreshold,clusterdef);
    elseif strcmp(posneg,'both') == 1
        [TempClustersP,UncorrectedThresholdCoordsP,clustersizesP,allcoordsP,allvalsP]=ClusterIdentifierSimple(vals,maskcoords,threshtype,thresholds(i,1),clustersizethreshold,clusterdef);
        [TempClustersN,UncorrectedThresholdCoordsN,clustersizesN,allcoordsN,allvalsN]=ClusterIdentifierSimple(vals*-1,maskcoords,threshtype,thresholds(i,1),clustersizethreshold,clusterdef);
        if iscell(TempClustersN)==0
            TempClustersN={TempClustersN};
        end;
        if iscell(TempClustersP)==0
            TempClustersP={TempClustersP};
        end;    
        if isempty(TempClustersN{1,1})==1
            TempClusters=TempClustersP;
        elseif isempty(TempClustersP{1,1})==1
            TempClusters=TempClustersP;
        else
            TempClusters=[TempClustersP;TempClustersN];
        end;        
        UncorrectedThresholdCoords=[UncorrectedThresholdCoordsP;UncorrectedThresholdCoordsN];
        clustersizes=[clustersizesP;clustersizesN];
        allcoords=[allcoordsP;allcoordsN];
        allvals=[allvalsP;allvalsN];
     end;
     
     Clusters.ThresholdedClusters{1,i}=TempClusters;
     Clusters.AllThresholdedCoords{1,i}=UncorrectedThresholdCoords;
     Clusters.MaxClusterSizes(1,i)=max(clustersizes);
     Clusters.AllClusterSizes{1,i}=clustersizes;
    
end;
end;



