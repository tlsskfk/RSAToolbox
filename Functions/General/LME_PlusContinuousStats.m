function [Results] = LME_PlusContinuousStats(DataTable,lmem_formula,whichstats)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%{'lme_Estimate','lme_SE','lme_tStat','lme_DF','lme_pValue','lme_Lower','lme_Upper','lme_adjRsquared','pearson_r','pearson_p','spearman_r','spearman_p','reg_beta','reg_se','reg_t','reg_t_pval','reg_sse','reg_dfe','reg_dfr','reg_ssr','reg_f','reg_f_pval','reg_adjrsquare'}

if nargin == 2 
    whichstats={'regression','pearson_r','spearman_r'};
end
try
    lmer_Model=fitglme(DataTable,lmem_formula);
catch
    for i = 1:length(DataTable.Properties.VariableNames)
        if isnumeric(DataTable.(DataTable.Properties.VariableNames{i}))
        DataTable.(DataTable.Properties.VariableNames{i})=double(DataTable.(DataTable.Properties.VariableNames{i}));
        end
    end 
    lmer_Model=fitglme(DataTable,lmem_formula);
end
yVarName=lmer_Model.ResponseName;
xVarNames=lmer_Model.CoefficientNames;
xVarNames(:,1)=[];
xVarNames=xVarNames(:);
Results=dataset2table(lmer_Model.Coefficients);
Results.Properties.VariableNames={'Name','lme_Estimate','lme_SE','lme_tStat','lme_DF','lme_pValue','lme_Lower','lme_Upper'};
Results.Name{1,1}=strrep(Results.Name{1,1},'(Intercept)','Intercept');
Results.Properties.RowNames=Results.Name;
Results.Name=[]; 
Results=[Results,array2table(repmat(lmer_Model.Rsquared.Adjusted,[size(Results,1),1]),'VariableNames',{'lme_adjRsquared'})];        
tempCorr=[];
if any(contains(whichstats,'pearson_r'))
    Y=table2array(DataTable(:,yVarName));
    Xs=table2array(DataTable(:,xVarNames));
    Xs=[ones(size(Xs,1),1),Xs];
    [tempCorr.pearson_r,tempCorr.pearson_p]=corr(Xs,Y,'rows','pairwise');
end    
if any(contains(whichstats,'spearman_r'))
    Y=table2array(DataTable(:,yVarName));
    Xs=table2array(DataTable(:,xVarNames));
    Xs=[ones(size(Xs,1),1),Xs];
    [tempCorr.spearman_r,tempCorr.spearman_p]=corr(Xs,Y,'type','spearman','rows','pairwise');
end   
if ~isempty(tempCorr)
    Results=[Results,struct2table(tempCorr)];
end

if any(contains(whichstats,'regression'))
    Y=table2array(DataTable(:,yVarName));
    Xs=table2array(DataTable(:,xVarNames));
    rstats=regstats(Y,Xs,'linear',{'adjrsquare','tstat','fstat'});
    rstats.tstat.t_pval=rstats.tstat.pval;
    rstats.tstat=rmfield(rstats.tstat,{'dfe','pval'});
    tstats=struct2table(rstats.tstat);
    rstats.fstat.f_pval=rstats.fstat.pval;
    rstats.fstat=rmfield(rstats.fstat,{'pval'});
    rstats.fstat.adjrsquare=rstats.adjrsquare;
    fstats=repmat(struct2table(rstats.fstat),[size(tstats,1),1]);
    reg_table=[tstats,fstats];
    reg_table.Properties.VariableNames={'reg_beta','reg_se','reg_t','reg_t_pval','reg_sse','reg_dfe','reg_dfr','reg_ssr','reg_f','reg_f_pval','reg_adjrsquare'};
    Results=[Results,reg_table];
end
