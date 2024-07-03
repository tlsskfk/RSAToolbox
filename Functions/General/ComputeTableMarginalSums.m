function [NewT,NewTTranspose,groupResults] = ComputeTableMarginalSums(T)
%Written by David Rothlein
    VarNames=T.Properties.VariableNames;
    RowNames=T.Properties.RowNames;
    NewVarNames=[{'AbsSum','NetSum','ContinuousLouvainGroupNum','p50LouvainGroupNum','p25LouvainGroupNum','p10LouvainGroupNum','p05LouvainGroupNum','Abs_ContinuousLouvainGroupNum','Abs_p50LouvainGroupNum','Abs_p25LouvainGroupNum','Abs_p10LouvainGroupNum','Abs_p05LouvainGroupNum'},VarNames];
    NewRowNames=[{'AbsSum';'NetSum';'ContinuousLouvainGroupNum';'p50LouvainGroupNum';'p25LouvainGroupNum';'p10LouvainGroupNum';'p05LouvainGroupNum';'Abs_ContinuousLouvainGroupNum';'Abs_p50LouvainGroupNum';'Abs_p25LouvainGroupNum';'Abs_p10LouvainGroupNum';'Abs_p05LouvainGroupNum'};RowNames];
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
    AbsCorrMatA=abs(CorrMatA);
    AbsCorrMatAT=abs(CorrMatAT);
    
    M_A(SelectIndAT)=community_louvain(CorrMatA,1,[],'negative_sym')';    
    M_AT(SelectIndA)=community_louvain(CorrMatAT,1,[],'negative_sym')';
    M_A_abs(SelectIndAT)=community_louvain(AbsCorrMatA)';    
    M_AT_abs(SelectIndA)=community_louvain(AbsCorrMatAT)';    
    M_A_p50(SelectIndAT)=community_louvain((single(CorrMatA < prctile(mat2uppertriuvectormat(CorrMatA),25))*-1 + single(CorrMatA > prctile(mat2uppertriuvectormat(CorrMatA),75))),1,[],'negative_sym')';
    M_AT_p50(SelectIndA)=community_louvain((single(CorrMatAT < prctile(mat2uppertriuvectormat(CorrMatAT),25))*-1 + single(CorrMatAT > prctile(mat2uppertriuvectormat(CorrMatAT),75))),1,[],'negative_sym')'; 
    M_A_p25(SelectIndAT)=community_louvain((single(CorrMatA < prctile(mat2uppertriuvectormat(CorrMatA),12.5))*-1 + single(CorrMatA > prctile(mat2uppertriuvectormat(CorrMatA),87.5))),1,[],'negative_sym')';
    M_AT_p25(SelectIndA)=community_louvain((single(CorrMatAT < prctile(mat2uppertriuvectormat(CorrMatAT),12.5))*-1 + single(CorrMatAT > prctile(mat2uppertriuvectormat(CorrMatAT),87.5))),1,[],'negative_sym')';
    M_A_p10(SelectIndAT)=community_louvain((single(CorrMatA < prctile(mat2uppertriuvectormat(CorrMatA),5))*-1 + single(CorrMatA > prctile(mat2uppertriuvectormat(CorrMatA),95))),1,[],'negative_sym')';
    M_AT_p10(SelectIndA)=community_louvain((single(CorrMatAT < prctile(mat2uppertriuvectormat(CorrMatAT),5))*-1 + single(CorrMatAT > prctile(mat2uppertriuvectormat(CorrMatAT),95))),1,[],'negative_sym')';  
    M_A_p5(SelectIndAT)=community_louvain((single(CorrMatA < prctile(mat2uppertriuvectormat(CorrMatA),2.5))*-1 + single(CorrMatA > prctile(mat2uppertriuvectormat(CorrMatA),97.5))),1,[],'negative_sym')';
    M_AT_p5(SelectIndA)=community_louvain((single(CorrMatAT < prctile(mat2uppertriuvectormat(CorrMatAT),2.5))*-1 + single(CorrMatAT > prctile(mat2uppertriuvectormat(CorrMatAT),97.5))),1,[],'negative_sym')';     
    M_A_abs_p50(SelectIndAT)=community_louvain(single(AbsCorrMatA > prctile(mat2uppertriuvectormat(AbsCorrMatA),50)))';
    M_AT_abs_p50(SelectIndA)=community_louvain(single(AbsCorrMatAT > prctile(mat2uppertriuvectormat(AbsCorrMatAT),50)))';
    M_A_abs_p25(SelectIndAT)=community_louvain(single(AbsCorrMatA > prctile(mat2uppertriuvectormat(AbsCorrMatA),75)))';
    M_AT_abs_p25(SelectIndA)=community_louvain(single(AbsCorrMatAT > prctile(mat2uppertriuvectormat(AbsCorrMatAT),75)))'; 
    M_A_abs_p10(SelectIndAT)=community_louvain(single(AbsCorrMatA > prctile(mat2uppertriuvectormat(AbsCorrMatA),90)))';
    M_AT_abs_p10(SelectIndA)=community_louvain(single(AbsCorrMatAT > prctile(mat2uppertriuvectormat(AbsCorrMatAT),90)))';  
    M_A_abs_p05(SelectIndAT)=community_louvain(single(AbsCorrMatA > prctile(mat2uppertriuvectormat(AbsCorrMatA),95)))';
    M_AT_abs_p05(SelectIndA)=community_louvain(single(AbsCorrMatAT > prctile(mat2uppertriuvectormat(AbsCorrMatAT),95)))';    
    
    add_A=zeros(length(SelectIndAT)-length(M_A),1);
    add_AT=zeros(length(SelectIndA)-length(M_AT),1);    
    AT=[AbsSumA,SumA,[M_A(:);add_A],[M_A_p50(:);add_A],[M_A_p25(:);add_A],[M_A_p10(:);add_A],[M_A_p5(:);add_A],[M_A_abs(:);add_A],[M_A_abs_p50(:);add_A],[M_A_abs_p25(:);add_A],[M_A_abs_p10(:);add_A],[M_A_abs_p05(:);add_A],AT];
    A=[AbsSumAT,SumAT,[M_AT(:);add_AT],[M_AT_p50(:);add_AT],[M_AT_p25(:);add_AT],[M_AT_p10(:);add_AT],[M_AT_p5(:);add_AT],[M_AT_abs(:);add_AT],[M_AT_abs_p50(:);add_AT],[M_AT_abs_p25(:);add_AT],[M_AT_abs_p10(:);add_AT],[M_AT_abs_p05(:);add_AT],A];
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
    groupResults.AbsCorrMatA=AbsCorrMatA;
    groupResults.AbsCorrMatAT=AbsCorrMatAT;
    groupResults.M_Inds.M_A=M_A(SelectIndAT)';
    groupResults.M_Inds.M_A_abs=M_A_abs(SelectIndAT)';
    groupResults.M_Inds.M_A_p50=M_A_p50(SelectIndAT)';
    groupResults.M_Inds.M_A_p25=M_A_p25(SelectIndAT)';
    groupResults.M_Inds.M_A_p10=M_A_p10(SelectIndAT)';
    groupResults.M_Inds.M_A_p05=M_A_p5(SelectIndAT)';
    groupResults.M_Inds.M_A_abs_p50=M_A_abs_p50(SelectIndAT)';
    groupResults.M_Inds.M_A_abs_p25=M_A_abs_p25(SelectIndAT)';
    groupResults.M_Inds.M_A_abs_p10=M_A_abs_p10(SelectIndAT)';
    groupResults.M_Inds.M_A_abs_p05=M_A_abs_p05(SelectIndAT)';
    
    groupResults.M_Inds.M_AT=M_AT(SelectIndA)';
    groupResults.M_Inds.M_AT_abs=M_AT_abs(SelectIndA)';
    groupResults.M_Inds.M_AT_p50=M_AT_p50(SelectIndA)';
    groupResults.M_Inds.M_AT_p25=M_AT_p25(SelectIndA)';
    groupResults.M_Inds.M_AT_p10=M_AT_p10(SelectIndA)';
    groupResults.M_Inds.M_AT_p05=M_AT_p5(SelectIndA)';
    groupResults.M_Inds.M_AT_abs_p50=M_AT_abs_p50(SelectIndA)';
    groupResults.M_Inds.M_AT_abs_p25=M_AT_abs_p25(SelectIndA)';
    groupResults.M_Inds.M_AT_abs_p10=M_AT_abs_p10(SelectIndA)';
    groupResults.M_Inds.M_AT_abs_p05=M_AT_abs_p05(SelectIndA)';    
    
    [groupResults.groupLabels.M_A,groupResults.groupTables.M_A] = groupInd2groupLabels(M_A,NewTTranspose,CorrMatA);
    [groupResults.groupLabels.M_A_abs,groupResults.groupTables.M_A_abs] = groupInd2groupLabels(M_A_abs,NewTTranspose,AbsCorrMatA);
    [groupResults.groupLabels.M_A_p50,groupResults.groupTables.M_A_p50] = groupInd2groupLabels(M_A_p50,NewTTranspose,CorrMatA);
    [groupResults.groupLabels.M_A_p25,groupResults.groupTables.M_A_p25] = groupInd2groupLabels(M_A_p25,NewTTranspose,CorrMatA);
    [groupResults.groupLabels.M_A_p10,groupResults.groupTables.M_A_p10] = groupInd2groupLabels(M_A_p10,NewTTranspose,CorrMatA);
    [groupResults.groupLabels.M_A_p05,groupResults.groupTables.M_A_p05] = groupInd2groupLabels(M_A_p5,NewTTranspose,CorrMatA);
    [groupResults.groupLabels.M_A_abs_p50,groupResults.groupTables.M_A_abs_p50] = groupInd2groupLabels(M_A_abs_p50,NewTTranspose,AbsCorrMatA);
    [groupResults.groupLabels.M_A_abs_p25,groupResults.groupTables.M_A_abs_p25] = groupInd2groupLabels(M_A_abs_p25,NewTTranspose,AbsCorrMatA);
    [groupResults.groupLabels.M_A_abs_p10,groupResults.groupTables.M_A_abs_p10] = groupInd2groupLabels(M_A_abs_p10,NewTTranspose,AbsCorrMatA);
    [groupResults.groupLabels.M_A_abs_p05,groupResults.groupTables.M_A_abs_p05] = groupInd2groupLabels(M_A_abs_p05,NewTTranspose,AbsCorrMatA);
 
    [groupResults.groupLabels.M_AT,groupResults.groupTables.M_AT] = groupInd2groupLabels(M_AT,NewT,CorrMatAT);
    [groupResults.groupLabels.M_AT_abs,groupResults.groupTables.M_AT_abs] = groupInd2groupLabels(M_AT_abs,NewT,AbsCorrMatAT);
    [groupResults.groupLabels.M_AT_p50,groupResults.groupTables.M_AT_p50] = groupInd2groupLabels(M_AT_p50,NewT,CorrMatAT);
    [groupResults.groupLabels.M_AT_p25,groupResults.groupTables.M_AT_p25] = groupInd2groupLabels(M_AT_p25,NewT,CorrMatAT);
    [groupResults.groupLabels.M_AT_p10,groupResults.groupTables.M_AT_p10] = groupInd2groupLabels(M_AT_p10,NewT,CorrMatAT);
    [groupResults.groupLabels.M_AT_p05,groupResults.groupTables.M_AT_p05] = groupInd2groupLabels(M_AT_p5,NewT,CorrMatAT);
    [groupResults.groupLabels.M_AT_abs_p50,groupResults.groupTables.M_AT_abs_p50] = groupInd2groupLabels(M_AT_abs_p50,NewT,AbsCorrMatAT);
    [groupResults.groupLabels.M_AT_abs_p25,groupResults.groupTables.M_AT_abs_p25] = groupInd2groupLabels(M_AT_abs_p25,NewT,AbsCorrMatAT);
    [groupResults.groupLabels.M_AT_abs_p10,groupResults.groupTables.M_AT_abs_p10] = groupInd2groupLabels(M_AT_abs_p10,NewT,AbsCorrMatAT);
    [groupResults.groupLabels.M_AT_abs_p05,groupResults.groupTables.M_AT_abs_p05] = groupInd2groupLabels(M_AT_abs_p05,NewT,AbsCorrMatAT);   

end

function [groupLabels,groupTables] = groupInd2groupLabels(groupInds,groupTable,corrMat)
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
        NetFeatSum=abs(nansum(tempFeatMat,1));
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
