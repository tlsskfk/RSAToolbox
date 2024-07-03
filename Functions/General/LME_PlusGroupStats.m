function [Results,PairwisePermTable,IndvPermTable] = LME_PlusGroupStats(DataTable,lmem_formula,groupLabels,PermDiffs,PermReps)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%{'lme_Estimate','lme_SE','lme_tStat','lme_DF','lme_pValue','lme_Lower','lme_Upper','lme_adjRsquared','pearson_r','pearson_p','spearman_r','spearman_p','reg_beta','reg_se','reg_t','reg_t_pval','reg_sse','reg_dfe','reg_dfr','reg_ssr','reg_f','reg_f_pval','reg_adjrsquare'}
numGroups=max(DataTable.GroupVar);
if nargin == 2 
    groupLabels=[];
    for i = 1:numGroups
        groupLabels{i,1}=num2str(i);
    end
    PermDiffs=1;
    PermReps=10000;
end
if nargin == 3 
    PermDiffs=1;
    PermReps=10000;
end
if nargin == 4
    PermReps=10000;
end
[~,groupPairLabels,~,~,PairNodes]=labels2uppertriuvectorlabels(groupLabels);
groupPairLabels=strrepCell(groupPairLabels,'_2_','vs');
DataTable.GroupVar=categorical(DataTable.GroupVar);
lmer_Model=fitglme(DataTable,lmem_formula);
yVarName=lmer_Model.ResponseName;
xVarNames=lmer_Model.CoefficientNames;
xVarNames(:,1)=[];
xVarNames=xVarNames(:);
Results=dataset2table(lmer_Model.Coefficients);
Results.Properties.VariableNames={'Name','lme_Estimate','lme_SE','lme_tStat','lme_DF','lme_pValue','lme_Lower','lme_Upper'};
Results.Name{1,1}=strrep(Results.Name{1,1},'(Intercept)','Intercept');
Results.Name{2,1}=strrep(Results.Name{2,1},'_2','');
Results.Properties.RowNames=Results.Name;
Results.Name=[]; 
Results=[Results,array2table(repmat(lmer_Model.Rsquared.Adjusted,[size(Results,1),1]),'VariableNames',{'lme_adjRsquared'})];        

if PermDiffs==1
    DataTable.GroupVar=single(DataTable.GroupVar);
    GroupData=cell(numGroups,1);
    for i = 1:numGroups
        GroupData{i,1}=DataTable.(yVarName)(DataTable.GroupVar==i,1);
    end    
    tableLabels={'Diff','PermP','PermCI','Tval','Tp','TCI'};
    [PermResults] = PairwiseDiffPermStats_noPar(GroupData,PermReps,0);
    PairwisePermTable=[];
    IndvPermTable=[];
    for i = 1:size(groupPairLabels,1)
        m=PairNodes(i,1);
        n=PairNodes(i,2);
        PairwisePermTable=[PairwisePermTable;{PermResults.PairwiseRealDiff(m,n),PermResults.PairwisePermPs(m,n),PermResults.PairwisePermCI{m,n},PermResults.PairwiseTVal(m,n),PermResults.PairwiseTPs(m,n),PermResults.PairwiseTCI{m,n}}];    
    end
    PairwisePermTable=cell2table(PairwisePermTable,'VariableNames',tableLabels,'RowNames',groupPairLabels);
    for i = 1:size(groupLabels,1)
        IndvPermTable=[IndvPermTable;{PermResults.IndvRealDiff(i,1),PermResults.IndvPermPs(i,1),PermResults.IndvPermCI{i,1},PermResults.IndvTVal(i,1),PermResults.IndvTPs(i,1),PermResults.IndvTCI{i,1}}]; 
    end    
    IndvPermTable=cell2table(IndvPermTable,'VariableNames',tableLabels,'RowNames',groupLabels);
end
end

