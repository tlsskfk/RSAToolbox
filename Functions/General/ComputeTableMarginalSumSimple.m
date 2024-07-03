function [NewT,NewTTranspose,groupResults] = ComputeTableMarginalSumSimple(T,AbsVal,Thresh,sortType)
%Written by David Rothlein
    if nargin == 1
        AbsVal=0;
        Thresh=0;
        sortType = 'abs';
    end    
    if nargin == 2
        Thresh = 0;
        sortType = 'abs';
    end
    if nargin == 3
        sortType = 'abs';
    end
    if Thresh > 0 && Thresh < 1
        Thresh = round(Thresh*100);
    end
    VarNames=T.Properties.VariableNames;
    RowNames=T.Properties.RowNames;
    if isempty(RowNames)
        RowNames=cell(height(T),1);
        for i = 1:height(T)
            RowNames{i,1}=['r',num2str(i)];
        end
    end
    NewVarNames=[{'AbsSum','NetSum','LouvainGroupNum'},VarNames];
    NewRowNames=[{'AbsSum';'NetSum';'LouvainGroupNum'};RowNames];
    A=table2array(T);
    AT=A';
    AbsSumA=nansum(abs(AT),2);
    AbsSumAT=nansum(abs(A),2);
    SelectIndAT=AbsSumA~=0;
    SelectIndA=AbsSumAT~=0;
    SumA=nansum(AT,2);
    SumAT=nansum(A,2);    
    CorrMatA=corr(A(:,SelectIndAT),'rows','pairwise');
    CorrMatAT=corr(AT(:,SelectIndA),'rows','pairwise');
    if AbsVal==1
        AbsCorrMatA=abs(CorrMatA);
        if Thresh ~= 0
            CorrMatA = community_louvain(single(CorrMatA > prctile(mat2uppertriuvectormat(AbsCorrMatA),100-Thresh)))';
        end
        M_A(SelectIndAT)=community_louvain(CorrMatA)'; 
    else
        if Thresh ~= 0
            CorrMatA = (single(CorrMatA < prctile(mat2uppertriuvectormat(CorrMatA),round(Thresh/2)))*-1 + single(CorrMatA > prctile(mat2uppertriuvectormat(CorrMatA),(100-round(Thresh/2)))));
        end
        M_A(SelectIndAT)=community_louvain(CorrMatA,1,[],'negative_sym')';    
    end
    if AbsVal==1
        AbsCorrMatAT=abs(CorrMatAT);
        if Thresh ~= 0
            CorrMatAT = community_louvain(single(CorrMatAT > prctile(mat2uppertriuvectormat(AbsCorrMatAT),100-Thresh)))';
        end
        M_AT(SelectIndA)=community_louvain(CorrMatAT)'; 
    else
        if Thresh ~= 0
            CorrMatAT = (single(CorrMatAT < prctile(mat2uppertriuvectormat(CorrMatAT),round(Thresh/2)))*-1 + single(CorrMatAT > prctile(mat2uppertriuvectormat(CorrMatAT),(100-round(Thresh/2)))));
        end
        M_AT(SelectIndA)=community_louvain(CorrMatAT,1,[],'negative_sym')';    
    end        
    add_A=zeros(length(SelectIndAT)-length(M_A),1);
    add_AT=zeros(length(SelectIndA)-length(M_AT),1);    
    AT=[AbsSumA,SumA,[M_A(:);add_A],AT];
    A=[AbsSumAT,SumAT,[M_AT(:);add_AT],A];
    NewT=array2table(A,'VariableNames',NewVarNames,'RowNames',RowNames);
    try
        NewTTranspose=array2table(AT,'VariableNames',NewRowNames,'RowNames',VarNames);   
    catch
        for n = 1:length(NewRowNames)
            if length(NewRowNames{n}) > 63
                NewRowNames{n}=NewRowNames{n}(1,1:63);
            end
        end
        NewTTranspose=array2table(AT,'VariableNames',NewRowNames,'RowNames',VarNames);  
    end
    
    groupResults.CorrMatA=CorrMatA;
    groupResults.CorrMatAT=CorrMatAT;
    groupResults.M_Inds.M_A=M_A(SelectIndAT)';
    groupResults.M_Inds.M_AT=M_AT(SelectIndA)';    
    
    [groupResults.groupLabels_M_A,groupResults.groupTables_M_A] = groupInd2groupLabels(M_A,NewTTranspose,CorrMatA,sortType);
    [groupResults.groupLabels_M_AT,groupResults.groupTables_M_AT] = groupInd2groupLabels(M_AT,NewT,CorrMatAT,sortType);
end

function [groupLabels,groupTables] = groupInd2groupLabels(groupInds,groupTable,corrMat,sortType)
    numGroups=max(groupInds);
    [~,modeSize]=mode(groupInds(:));
    groupLabels=cell(modeSize,numGroups);
    groupTables=cell(1,numGroups);
    cleanInds=groupInds(groupInds~=0);
    LabelNames=groupTable.Properties.RowNames;
    groupTable(:,3:12)=[];
    sumTable=groupTable(:,[1,2]);
    featTable=groupTable(:,3:end);
    featMat=table2array(featTable);
    for i = 1:numGroups        
        tempMat=corrMat;
        temp=LabelNames(groupInds==i);
        groupLabels(1:length(temp),i)= temp(:);
        tempFeatMat=featMat(groupInds==i,:);
        FeatSelect=nansum(abs(tempFeatMat),1)>0;
        tempFeatMat=tempFeatMat(:,FeatSelect);
        tempCorrMat=corr(tempFeatMat');
        signInd=single(sign(tempCorrMat(:,1)));
        if sum(signInd,1) < 0
            signInd=signInd*-1;
        end
        tempFeatMat=tempFeatMat.*repmat(signInd,[1,size(tempFeatMat,2)]);
        tempFeatTable=featTable(groupInds==i,FeatSelect);
        tempFeatTable=array2table(tempFeatMat,'VariableNames',tempFeatTable.Properties.VariableNames,'RowNames',tempFeatTable.Properties.RowNames);
        if strcmpi(sortType,'abs')
            NetFeatSum=abs(nansum(tempFeatMat,1));
        else
            NetFeatSum=nansum(tempFeatMat,1);
        end
        [~,ord]=sort(NetFeatSum,'descend');
        tempFeatTable=tempFeatTable(:,ord);
        groupTables{1,i}=[sumTable(groupInds==i,:),tempFeatTable]; 
        groupTables{1,i}.NetSum=nanmean(tempFeatMat,2);
        groupTables{1,i}.AbsSum=nanmean(abs(tempFeatMat),2);
        groupTables{1,i}.VarSign=signInd;
        groupTables{1,i} = movevars(groupTables{1,i}, 'VarSign', 'After', 'NetSum');
        groupCorrMat=corrMat(cleanInds==i,cleanInds==i);
        groupCorrMat(eye(length(groupCorrMat))==1)=nan;
        groupCorrMeans=nanmean(groupCorrMat,2);  
        groupTables{1,i}=[array2table(groupCorrMeans,'VariableNames',{'MeanWithinGroupCorr'}),groupTables{1,i}];
        groupTables{1,i}=[groupTables{1,i};array2table([nanmean(abs(table2array(groupTables{1,i})),1);nanmean(table2array(groupTables{1,i}),1)],'VariableNames',groupTables{1,i}.Properties.VariableNames,'RowNames',{'AbsFeatMean','FeatMean'})];        
        groupTables{2,i}=nanmean(mat2uppertriuvectormat(corrMat(cleanInds==i,cleanInds==i)));
        tempMat(cleanInds==i,cleanInds==i)=nan;
        groupTables{3,i}=nanmean(mat2uppertriuvectormat(tempMat));
        groupTables{4,i}=groupTables{2,i}-groupTables{3,i}; 
        groupTables{1,i}.MeanWithinGroupCorr([end-1,end],1)=nan(2,1);
        groupTables{1,i}.AbsSum([end-1,end],1)=nan(2,1);
        groupTables{1,i}.NetSum([end-1,end],1)=nan(2,1);
        groupTables{1,i}.VarSign([end-1,end],1)=nan(2,1);
    end
end
