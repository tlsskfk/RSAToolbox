function [ Reg_RTz_byItem_Grp,Reg_CEz_byItem_Grp,Reg_RTz_byItem_SS,Reg_RTz_byTrial_SS,Reg_CEz_byTrial_SS,IndvDiff_Reg_CEz_byTrial,IndvDiff_Reg_RTz_byTrial,IndvDiff_Corr_CEz_byTrial,IndvDiff_Corr_RTz_byTrial,RTz_Index_b,RTz_Index_r,RTz_Index_t, CEz_Index_b,CEz_Index_t,CEz_Index_r] = GradCPT_RegressAnalysis( RT_byItem_bySS, RT_byTrial_bySS, CE_byItem_Grp, CE_byTrial_bySS,ind_var_city,ind_var_mount,RT_Items,CE_Items,stats,linear)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Reg_RTz_byItem_Grp=[];
Reg_CEz_byItem_Grp=[];
Reg_RTz_byItem_SS=[];
Reg_RTz_byTrial_SS=[];
Reg_CEz_byTrial_SS=[];
IndvDiff_Reg_CEz_byTrial=[];
IndvDiff_Reg_RTz_byTrial=[];
IndvDiff_Corr_CEz_byTrial=[];
IndvDiff_Corr_RTz_byTrial=[];
RTz_Index_b=[];
RTz_Index_r=[];
RTz_Index_t=[];
CEz_Index_b=[];
CEz_Index_t=[];
CEz_Index_r=[];
Reg_CEz_byTrial_SS=[];
IndvDiff_Reg_CEz_byTrial=[];
IndvDiff_Corr_CEz_byTrial=[];
CEz_Index_b=[];
CEz_Index_t=[];
CEz_Index_r=[];
for i =1:size(RT_byItem_bySS,2)
    RTz_byItem_bySS(:,i)=nan_zscore(RT_byItem_bySS(:,i));
end
RT_byItem_Grp=nanmean(RTz_byItem_bySS,2);
RTz_byItem_Grp=nan_zscore(RT_byItem_Grp);

for i = 1:size(RT_byTrial_bySS,1)
    RTz_byTrial_bySS{i,1}=nan_zscore(RT_byTrial_bySS{i,1});
    CEz_byTrial_bySS{i,1}=nan_zscore(CE_byTrial_bySS{i,1});
end
    
CEz_byItem_Grp=nan_zscore(CE_byItem_Grp);


for i = 1:size(ind_var_city,2)
    ind_var_city_z(:,i)=zscore(ind_var_city(:,i));
    ind_var_mount_z(:,i)=zscore(ind_var_mount(:,i));
end

if sum(single(isnan(Reg_CEz_byItem_Grp)))==0
    Reg_RTz_byItem_Grp=regstats(RTz_byItem_Grp,ind_var_city_z,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
end
% if sum(single(isnan(Reg_CEz_byItem_Grp)))==0
% %    Reg_CEz_byItem_Grp=regstats(CEz_byItem_Grp,ind_var_mount_z,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
% end
for i = 1:size(RT_byTrial_bySS,1)
   % i/size(RT_byTrial_bySS,1)
    city_ind_var_ByTrial_z=[];
    mount_ind_var_ByTrial_z=[];
    RT_Item=RT_Items{i,1};
    CE_Item=CE_Items{i,1};
    if strcmp(linear,'linear')==1
        RTz_linearTrend=zscore([1:size(RT_Item,1)])';
        CEz_linearTrend=zscore([1:size(CE_Item,1)])';
    end
    
    city_ind_var_ByTrial=ind_var_city(RT_Item,:);
    mount_ind_var_ByTrial=ind_var_mount(CE_Item,:);
    
    for j = 1:size(ind_var_city,2)
        city_ind_var_ByTrial_z(:,j)=nan_zscore(city_ind_var_ByTrial(:,j));
        mount_ind_var_ByTrial_z(:,j)=nan_zscore(mount_ind_var_ByTrial(:,j));
    end
    
    
    tReg=regstats(RTz_byItem_bySS(:,i),ind_var_city_z,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
    tCorr=nancorr(RTz_byItem_bySS(:,i),ind_var_city_z);
    Reg_RTz_byItem_SS.Corr(i,:)=tCorr;
    Reg_RTz_byItem_SS.AdjR(i,:)=tReg.adjrsquare;
    Reg_RTz_byItem_SS.R(i,:)=tReg.rsquare;
    Reg_RTz_byItem_SS.f(i,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
    Reg_RTz_byItem_SS.t(i,:)=tReg.tstat.t(2:end);
    Reg_RTz_byItem_SS.b(i,:)=tReg.tstat.beta(2:end);
    Reg_RTz_byItem_SS.se(i,:)=tReg.tstat.se(2:end);
    Reg_RTz_byItem_SS.pval(i,:)=tReg.tstat.pval(2:end);
    Reg_RTz_byItem_SS.df(i,:)=tReg.tstat.dfe;
    
    if strcmp(linear,'linear')==1
        tReg=regstats(RTz_byTrial_bySS{i,1},[city_ind_var_ByTrial_z,RTz_linearTrend],'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
        tCorr=nancorr(RTz_byTrial_bySS{i,1},[city_ind_var_ByTrial_z,RTz_linearTrend]);
    else
        tReg=regstats(RTz_byTrial_bySS{i,1},[city_ind_var_ByTrial_z],'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
        tCorr=nancorr(RTz_byTrial_bySS{i,1},[city_ind_var_ByTrial_z]);  
    end
    
    Reg_RTz_byTrial_SS.Corr(i,:)=tCorr;
    Reg_RTz_byTrial_SS.AdjR(i,:)=tReg.adjrsquare;
    Reg_RTz_byTrial_SS.R(i,:)=tReg.rsquare;
    Reg_RTz_byTrial_SS.f(i,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
    Reg_RTz_byTrial_SS.t(i,:)=tReg.tstat.t(2:end);
    Reg_RTz_byTrial_SS.b(i,:)=tReg.tstat.beta(2:end);
    Reg_RTz_byTrial_SS.se(i,:)=tReg.tstat.se(2:end);
    Reg_RTz_byTrial_SS.pval(i,:)=tReg.tstat.pval(2:end);
    Reg_RTz_byTrial_SS.df(i,:)=tReg.tstat.dfe;
    
    if strcmp(linear,'linear')==1
        numX=size(mount_ind_var_ByTrial_z,2)+1;
    else
        numX=size(mount_ind_var_ByTrial_z,2);
    end
    
    if length(CEz_byTrial_bySS{i,1})<=numX || isnan(CEz_byTrial_bySS{i,1}(1,1))==1

        Reg_CEz_byTrial_SS.Corr(i,:)=nan(1,numX);
        Reg_CEz_byTrial_SS.AdjR(i,:)=nan;
        Reg_CEz_byTrial_SS.R(i,:)=nan;
        Reg_CEz_byTrial_SS.f(i,:)=nan(1,4);
        Reg_CEz_byTrial_SS.t(i,:)=nan(1,numX);
        Reg_CEz_byTrial_SS.b(i,:)=nan(1,numX);
        Reg_CEz_byTrial_SS.se(i,:)=nan(1,numX);
        Reg_CEz_byTrial_SS.pval(i,:)=nan(1,numX);
        Reg_CEz_byTrial_SS.df(i,:)=nan;
    else   
        if strcmp(linear,'linear')==1
            tReg=regstats(CEz_byTrial_bySS{i,1},[mount_ind_var_ByTrial_z,CEz_linearTrend],'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
            tCorr=nancorr(CEz_byTrial_bySS{i,1},[mount_ind_var_ByTrial_z,CEz_linearTrend]);
        else
            tReg=regstats(CEz_byTrial_bySS{i,1},[mount_ind_var_ByTrial_z],'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
            tCorr=nancorr(CEz_byTrial_bySS{i,1},[mount_ind_var_ByTrial_z]);
        end
                
        Reg_CEz_byTrial_SS.Corr(i,:)=tCorr;
        Reg_CEz_byTrial_SS.AdjR(i,:)=tReg.adjrsquare;
        Reg_CEz_byTrial_SS.R(i,:)=tReg.rsquare;
        Reg_CEz_byTrial_SS.f(i,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
        Reg_CEz_byTrial_SS.t(i,:)=tReg.tstat.t(2:end);
        Reg_CEz_byTrial_SS.b(i,:)=tReg.tstat.beta(2:end);
        Reg_CEz_byTrial_SS.se(i,:)=tReg.tstat.se(2:end);
        Reg_CEz_byTrial_SS.pval(i,:)=tReg.tstat.pval(2:end);
        Reg_CEz_byTrial_SS.df(i,:)=tReg.tstat.dfe;
    end
end

RTz_Index_b=nanmean(Reg_RTz_byTrial_SS.b(:,[1:2]),2);
RTz_Index_t=nanmean(Reg_RTz_byTrial_SS.t(:,[1:2]),2);
RTz_Index_r=nanmean(Reg_RTz_byTrial_SS.Corr(:,[1:2]),2);

CEz_Index_b=nanmean(Reg_CEz_byTrial_SS.b(:,[1:2]),2);
CEz_Index_t=nanmean(Reg_CEz_byTrial_SS.t(:,[1:2]),2);
CEz_Index_r=nanmean(Reg_CEz_byTrial_SS.Corr(:,[1:2]),2);
if size(stats,1)>1
    for j=1:size(stats,2)
        stat=stats(:,j);
        nUse=sum(sum(isnan(Reg_CEz_byTrial_SS.t)==0,2)==numX);
        if nUse <=numX
            IndvDiff_Reg_CEz_byTrial=nan;
            IndvDiff_Corr_CEz_byTrial=nan;
        else
            tReg=regstats(stat,Reg_CEz_byTrial_SS.t,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
            tCorr=nancorr(stat,[Reg_CEz_byTrial_SS.t,Reg_CEz_byTrial_SS.R]);
            [r_ind_r, p_ind_r]=nancorr(stat,CEz_Index_r);
            [r_ind_t, p_ind_t]=nancorr(stat,CEz_Index_t);
            [r_ind_b, p_ind_b]=nancorr(stat,CEz_Index_b);
            IndvDiff_Reg_CEz_byTrial.Corr(j,:)=tCorr;
            IndvDiff_Reg_CEz_byTrial.AdjR(j,:)=tReg.adjrsquare;
            IndvDiff_Reg_CEz_byTrial.R(j,:)=tReg.rsquare;
            IndvDiff_Reg_CEz_byTrial.f(j,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
            IndvDiff_Reg_CEz_byTrial.t(j,:)=tReg.tstat.t(2:end);
            IndvDiff_Reg_CEz_byTrial.b(j,:)=tReg.tstat.beta(2:end);
            IndvDiff_Reg_CEz_byTrial.se(j,:)=tReg.tstat.se(2:end);
            IndvDiff_Reg_CEz_byTrial.pval(j,:)=tReg.tstat.pval(2:end);
            IndvDiff_Reg_CEz_byTrial.df(j,:)=tReg.tstat.dfe;
            IndvDiff_Reg_CEz_byTrial.Ind_r.r(j,:)=r_ind_r;
            IndvDiff_Reg_CEz_byTrial.Ind_r.p(j,:)=p_ind_r;
            IndvDiff_Reg_CEz_byTrial.Ind_t.r(j,:)=r_ind_t;
            IndvDiff_Reg_CEz_byTrial.Ind_t.p(j,:)=p_ind_t;
            IndvDiff_Reg_CEz_byTrial.Ind_b.r(j,:)=r_ind_b;
            IndvDiff_Reg_CEz_byTrial.Ind_b.p(j,:)=p_ind_b;
            tReg=regstats(stat,Reg_CEz_byTrial_SS.Corr,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
            tCorr=nancorr(stat,[Reg_CEz_byTrial_SS.Corr,Reg_CEz_byTrial_SS.R]);
            IndvDiff_Corr_CEz_byTrial.Corr(j,:)=tCorr;
            IndvDiff_Corr_CEz_byTrial.AdjR(j,:)=tReg.adjrsquare;
            IndvDiff_Corr_CEz_byTrial.R(j,:)=tReg.rsquare;
            IndvDiff_Corr_CEz_byTrial.f(j,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
            IndvDiff_Corr_CEz_byTrial.t(j,:)=tReg.tstat.t(2:end);
            IndvDiff_Corr_CEz_byTrial.b(j,:)=tReg.tstat.beta(2:end);
            IndvDiff_Corr_CEz_byTrial.se(j,:)=tReg.tstat.se(2:end);
            IndvDiff_Corr_CEz_byTrial.pval(j,:)=tReg.tstat.pval(2:end);
            IndvDiff_Corr_CEz_byTrial.df(j,:)=tReg.tstat.dfe;
        end
        stat=nan_zscore(stats(:,j));
        tReg=regstats(stat,Reg_RTz_byTrial_SS.t,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
        tCorr=nancorr(stat,[Reg_RTz_byTrial_SS.t,Reg_RTz_byTrial_SS.R]);
        [r_ind_r, p_ind_r]=nancorr(stat,RTz_Index_r);
        [r_ind_t, p_ind_t]=nancorr(stat,RTz_Index_t);
        [r_ind_b, p_ind_b]=nancorr(stat,RTz_Index_b);
        IndvDiff_Reg_RTz_byTrial.Corr(j,:)=tCorr;
        IndvDiff_Reg_RTz_byTrial.AdjR(j,:)=tReg.adjrsquare;
        IndvDiff_Reg_RTz_byTrial.R(j,:)=tReg.rsquare;
        IndvDiff_Reg_RTz_byTrial.f(j,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
        IndvDiff_Reg_RTz_byTrial.t(j,:)=tReg.tstat.t(2:end);
        IndvDiff_Reg_RTz_byTrial.b(j,:)=tReg.tstat.beta(2:end);
        IndvDiff_Reg_RTz_byTrial.se(j,:)=tReg.tstat.se(2:end);
        IndvDiff_Reg_RTz_byTrial.pval(j,:)=tReg.tstat.pval(2:end);
        IndvDiff_Reg_RTz_byTrial.df(j,:)=tReg.tstat.dfe;        
        IndvDiff_Reg_RTz_byTrial.Ind_r.r(j,:)=r_ind_r;
        IndvDiff_Reg_RTz_byTrial.Ind_r.p(j,:)=p_ind_r;
        IndvDiff_Reg_RTz_byTrial.Ind_t.r(j,:)=r_ind_t;
        IndvDiff_Reg_RTz_byTrial.Ind_t.p(j,:)=p_ind_t;
        IndvDiff_Reg_RTz_byTrial.Ind_b.r(j,:)=r_ind_b;
        IndvDiff_Reg_RTz_byTrial.Ind_b.p(j,:)=p_ind_b;
        tReg=regstats(stat,Reg_RTz_byTrial_SS.Corr,'linear',{'tstat','rsquare','fstat','adjrsquare','standres'});
        tCorr=nancorr(stat,[Reg_RTz_byTrial_SS.Corr,Reg_RTz_byTrial_SS.R]);
        IndvDiff_Corr_RTz_byTrial.Corr(j,:)=tCorr;
        IndvDiff_Corr_RTz_byTrial.AdjR(j,:)=tReg.adjrsquare;
        IndvDiff_Corr_RTz_byTrial.R(j,:)=tReg.rsquare;
        IndvDiff_Corr_RTz_byTrial.f(j,:)=[tReg.fstat.dfe,tReg.fstat.dfr,tReg.fstat.f,tReg.fstat.pval];
        IndvDiff_Corr_RTz_byTrial.t(j,:)=tReg.tstat.t(2:end);
        IndvDiff_Corr_RTz_byTrial.b(j,:)=tReg.tstat.beta(2:end); 
        IndvDiff_Corr_RTz_byTrial.se(j,:)=tReg.tstat.se(2:end);
        IndvDiff_Corr_RTz_byTrial.pval(j,:)=tReg.tstat.pval(2:end);
        IndvDiff_Corr_RTz_byTrial.df(j,:)=tReg.tstat.dfe; 
    end
end


