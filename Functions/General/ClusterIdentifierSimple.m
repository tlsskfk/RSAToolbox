function [Clusters,UncorrectedThresholdCoords,clustersizes,allcoords,allvals]=ClusterIdentifierSimple(vals,maskcoords,threshtype,threshold,clustersizethreshold,clusterdef) %'side';'line';'point
%Written by David Rothlein
allcoords=[];
allvals=[];
maskcoords=single(maskcoords);
clustersearch1=maskcoords;
vals(isnan(vals)==1,:)=0;
vals(isinf(vals)==1,:)=0;

if strcmp(threshtype,'rank')==1
    [~,rankindex]=sort(vals);
    if threshold<1
        threshold = round(threshold*length(vals));
    end
    ThreshFilter=rankindex(1:(size(rankindex,1)-threshold));
else
    ThreshFilter=vals<threshold;
end

clustersearch1(ThreshFilter,:)=[];
vals(ThreshFilter,:)=[];
UncorrectedThresholdCoords=clustersearch1;
clustersearch2=[clustersearch1,zeros(size(clustersearch1,1),1)];

if size(clustersearch1,1) == 0
    UncorrectedThresholdCoords=[];
    Clusters=[];
    clustersizes=0;
else
    endcount=0;
    clustercount=1;
    
    while endcount<1
        clustersearch1=UncorrectedThresholdCoords;
        tempfilter = clustersearch2(:,4)>0;
        clustersearch1(tempfilter,:)=[];
        
        [ind,~]=coord2cluster(clustersearch1(1,:),UncorrectedThresholdCoords,clusterdef);
        clustersearch2(:,4)=clustersearch2(:,4)+(ind*clustercount);
        clustercount=clustercount+1;
        clustersearch1=UncorrectedThresholdCoords;
        tempfilter = clustersearch2(:,4)>0;
        clustersearch1(tempfilter,:)=[];        
        
        if size(clustersearch1,1)==0
            endcount = 1;
        end
    end
    
    numofclusters=max(clustersearch2(:,4));
    clustersizes=zeros(numofclusters,1);
    Clusters=cell(1,3);
    count=1;
    allcoords=[];
    allvals=[];
    for i=1:numofclusters
        
        clustfilter=clustersearch2(:,4)==i;
        clustersizes(i,1)=sum(clustfilter(:));
        if clustersizes(i,1)>=clustersizethreshold
            Clusters{count,1}=UncorrectedThresholdCoords(clustfilter,:);
            Clusters{count,2}=vals(clustfilter,:);
            Clusters{count,3}=clustersizes(i,1); 

            count=count+1;           
        end
        
    end
    for i=1:size(Clusters,1)
        allcoords=[allcoords;Clusters{i,1}];
        allvals=[allvals;Clusters{i,2}];
    end
end
end


function [ind,clustercoords]=coord2cluster(seedcoords,coords,clusterdef) %'side';'line';'point'

ind=zeros(size(coords,1),1);
used=ind+1;

tempsearch=[abs(coords(:,1)-seedcoords(1,1)),abs(coords(:,2)-seedcoords(1,2)),abs(coords(:,3)-seedcoords(1,3))];

if strcmpi(clusterdef,'side')==1 || strcmpi(clusterdef,'face')==1
    tempfilter=sum(tempsearch<2,2)==3;
    tempfilter2=sum(tempsearch==0,2)>1;
    tempfilter=tempfilter.*tempfilter2;
    tempfilter=tempfilter>0;
end

if strcmpi(clusterdef,'line')==1 || strcmpi(clusterdef,'edge')==1
    tempfilter=sum(tempsearch<2,2)==3;
    tempfilter2=sum(tempsearch==0,2)>0;
    tempfilter=tempfilter.*tempfilter2;
    tempfilter=tempfilter>0;
end

if strcmpi(clusterdef,'point')==1 || strcmpi(clusterdef,'vertex')==1 || strcmpi(clusterdef,'node')==1
    tempfilter=sum(tempsearch<2,2)==3;
end

ind(tempfilter,1)=1;
used(tempfilter,1)=0;
oldsum=sum(ind(:));
count=1;

if oldsum == 1
    newsum=0;
    clustercoords=seedcoords;
else    
    newsum=oldsum+1;
    seedcoords=coords((ind>0),:);
end

while oldsum>0
    oldsum=size(seedcoords,1);    
    for i = 1:oldsum
        tempsearch=[abs(coords(:,1)-seedcoords(i,1)),abs(coords(:,2)-seedcoords(i,2)),abs(coords(:,3)-seedcoords(i,3))]; 
        
        if strcmpi(clusterdef,'side')==1 || strcmpi(clusterdef,'face')==1
            tempfilter=sum(tempsearch<2,2)==3;
            tempfilter2=sum(tempsearch==0,2)>1;
            tempfilter=tempfilter.*tempfilter2;
            tempfilter=tempfilter>0;
        end

        if strcmpi(clusterdef,'line')==1 || strcmpi(clusterdef,'edge')==1
            tempfilter=sum(tempsearch<2,2)==3;
            tempfilter2=sum(tempsearch==0,2)>0;
            tempfilter=tempfilter.*tempfilter2;
            tempfilter=tempfilter>0;
        end

        if strcmpi(clusterdef,'point')==1 || strcmpi(clusterdef,'vertex')==1 || strcmpi(clusterdef,'node')==1
            tempfilter=sum(tempsearch<2,2)==3;
        end
        tempfilter=tempfilter>0&ind==0;
        %tempfilter=(sum(tempsearch<2,2)==3&ind==0);        
        
        
        ind(tempfilter,1)=count+1;
    end
    count=count+1;
    seedcoords=coords((ind==(count)),:);
    
    oldsum=size(seedcoords,1);
end
ind=ind>0;
clustercoords=coords((ind>0),:);
end