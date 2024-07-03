function [IndvDiffTables,bonfMat] = Feature2IndvDiffTables(FeatTables,BehTable,varargin)
%Written by David Rothlein
[StatType] = VariableSetter('StatType','pearson',varargin);
[Zfeat] = VariableSetter('Zfeat',0,varargin);
[centerfeat] = VariableSetter('centerfeat',0,varargin);
if iscell(FeatTables)
    numParcels=size(FeatTables,1);
    numFeatTypes=size(FeatTables,2);
else
    numParcels=1;
    numFeatTypes=1;
    FeatTables={FeatTables};
end
BehVarNames=BehTable.Properties.VariableNames;
BehTable=table2array(BehTable);
bonfMat=zeros(numParcels,numFeatTypes);
count=1;
for FeatNum = 1:numFeatTypes
    for parcelNum=1:numParcels
        disp([num2str(round((count/(numFeatTypes*numParcels))*100)),'% complete...'])
        count=count+1;
        FeatNames = FeatTables{parcelNum,FeatNum}.Properties.VariableNames;
        numFeature=length(FeatNames);
        FeatTable=table2array(FeatTables{parcelNum,FeatNum});
        if Zfeat == 1
            FeatTable=zscore(FeatTable')';
        end   
        if centerfeat == 1
            FeatTable=FeatTable-repmat(nanmean(FeatTable,2),[1,size(FeatTable,2)]);
        end          
        if strcmpi(StatType,'pearson')
            [tempStat,p]=corr(BehTable,FeatTable,'rows','pairwise');
        elseif strcmpi(StatType,'spearman')
            [tempStat,p]=corr(BehTable,FeatTable,'rows','pairwise','type','spearman');
        end
        IndvDiffTables.p{parcelNum,FeatNum}=array2table(SigSeg2(p,[0.1,0.05,0.01,0.005,0.001,0.0005,0.0001],'pos').*sign(tempStat),'VariableNames',FeatNames,'RowNames',BehVarNames);
        tempBonf = SigSeg2(p,0.05/numFeature,'pos');
        yRemove=sum(tempBonf,1)==0;
        xRemove=sum(tempBonf,2)==0;
        BonfSum=sum(tempBonf(:));
        tempBonf= array2table(tempBonf.*sign(tempStat),'VariableNames',FeatNames,'RowNames',BehVarNames);
        tempBonf(xRemove,:)=[];
        tempBonf(:,yRemove)=[];
        IndvDiffTables.bonf{parcelNum,FeatNum}=tempBonf;
        bonfMat(parcelNum,FeatNum)=BonfSum;
        IndvDiffTables.stat{parcelNum,FeatNum}=array2table(tempStat,'VariableNames',FeatNames,'RowNames',BehVarNames);
    end
end


end