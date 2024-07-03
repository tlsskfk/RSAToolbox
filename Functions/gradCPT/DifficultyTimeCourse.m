function [ DifficultyTC, MountainTC, CityTC, CleanVTC, Mountain_r, City_r,CityTCSplit,splitlabels,TC_summary] = DifficultyTimeCourse( CityTrials, MountainTrials, ind_var_city, ind_var_mount, CityItems, MountainItems, CityRTs, MountainCE, VTC_byTrial, CityRT_byItem )
DifficultyTC=[];
MountainTC=[];
CleanVTC=[];
Mountain_r=[];
CityTC=[];
MeanRTz=[-0.929401720966902;-0.517149778760764;-0.939049608822749;0.248996769957092;-0.701862128783757;2.03550475979877;0.917389758286787;-0.269194463956843;0.902919082210006;-0.748152668961644]; %computedFrom TMB dataset
splitlabels{1,1}='TargetDensity';
splitlabels{1,2}='Similarity';
splitlabels{1,3}='Exemplar';
splitlabels{1,4}='ExemplarMean';
splitlabels{1,5}='Category';
splitlabels{1,6}='ExemplarPix';
splitlabels{1,7}='ExemplarGist';
splitlabels{1,8}='RTz_TC_TMB';
splitlabels{1,9}='RTz_TC';
splitlabels{1,10}='RTz_TC_TMBxVTC';
splitlabels{1,11}='TargetDensityxVTC';
splitlabels{1,12}='RTs_demeaned_bySS';
splitlabels{1,13}='RTs_demeaned_Leave1OutGrp';
numofsubj=size(CityTrials,1);
ItemRTs_demeaned=CityRT_byItem-repmat(mean(CityRT_byItem),size(CityRT_byItem,1),1);
for i = 1:numofsubj
    TempCityTrials=CityTrials{i,1};
    TempCityTrials(end)=[];
    TempMountainTrials=MountainTrials{i,1};
    TempMountainTrials(end)=[];
    TempCityItems=CityItems{i,1};
    TempMountainItems=MountainItems{i,1};
    TempCityRTs=CityRTs{i,1};
    TempMountainCEs=MountainCE{i,1};
    TempVTC_byTrial=VTC_byTrial{i,1};
    Temp_ItemRTs_demeaned_SS=ItemRTs_demeaned(:,i);
    Temp_ItemRTs_demeaned_Grp=ItemRTs_demeaned;
    Temp_ItemRTs_demeaned_Grp(:,i)=[];
    Temp_ItemRTs_demeaned_Grp=nanmean(Temp_ItemRTs_demeaned_Grp,2);
    
    [ TempTargetDensityTC ] = TargetConcentrationTC( TempMountainTrials,0,14,1.5 ); %TargetConcentrationTC( TargetTrials,gauss_mu,gauss_sigma,z_ceiling )
    
%     CityTCSplit{i,1}=zscore(TempTargetDensityTC)';
%     CityTCSplit{i,2}=0.3141*zscore(ind_var_city(TempCityItems,1)) + 0.3333*zscore(ind_var_city(TempCityItems,2)) + 0.1640*zscore(ind_var_city(TempCityItems,3)) + -0.1214*zscore(ind_var_city(TempCityItems,4));
%     CityTCSplit{i,3}=0.3586*zscore(ind_var_city(TempCityItems,1)) + 0.1987*zscore(ind_var_city(TempCityItems,2));
%     CityTCSplit{i,4}=-0.0673*zscore(ind_var_city(TempCityItems,3)) + -0.2110*zscore(ind_var_city(TempCityItems,4));
%     CityTCSplit{i,5}=0.3578*zscore(ind_var_city(TempCityItems,1)) + -0.0944*zscore(ind_var_city(TempCityItems,3));
%     CityTCSplit{i,6}=0.3920*zscore(ind_var_city(TempCityItems,2)) + 0.2750*zscore(ind_var_city(TempCityItems,4));
%     CityTCSplit{i,7}=zscore(ind_var_city(TempCityItems,1));
%     CityTCSplit{i,8}=zscore(ind_var_city(TempCityItems,2));
%     CityTCSplit{i,9}=zscore(ind_var_city(TempCityItems,3));
%     CityTCSplit{i,10}=zscore(ind_var_city(TempCityItems,4));
%     CityTCSplit{i,12}=MeanRTz(TempCityItems,1);
    
    CityTCSplit{i,1}=zscore(TempTargetDensityTC)';
    CityTCSplit{i,2}=0.3402*zscore(ind_var_city(TempCityItems,1)) + 0.2282*zscore(ind_var_city(TempCityItems,2)) + 0.0699*zscore(ind_var_city(TempCityItems,3));
    CityTCSplit{i,3}=0.3586*zscore(ind_var_city(TempCityItems,1)) + 0.1987*zscore(ind_var_city(TempCityItems,2));
    CityTCSplit{i,4}=(zscore(ind_var_city(TempCityItems,1)) + zscore(ind_var_city(TempCityItems,2)))/2;
    CityTCSplit{i,5}=zscore(ind_var_city(TempCityItems,3));
    CityTCSplit{i,6}=zscore(ind_var_city(TempCityItems,1));
    CityTCSplit{i,7}=zscore(ind_var_city(TempCityItems,2));
    CityTCSplit{i,8}=MeanRTz(TempCityItems,1);    
    CityTCSplit{i,10}=CityTCSplit{i,8}.*TempVTC_byTrial(TempCityTrials);
    CityTCSplit{i,11}=CityTCSplit{i,1}.*TempVTC_byTrial(TempCityTrials);
    CityTCSplit{i,12}=Temp_ItemRTs_demeaned_SS(TempCityItems,1);  
    CityTCSplit{i,13}=Temp_ItemRTs_demeaned_Grp(TempCityItems,1);  
    TempCityRTs(isnan(TempCityRTs)==1)=nanmean(TempCityRTs);   
    CityTCSplit{i,9}=TempCityRTs;
    for j = 1:size(CityTCSplit,2)
        City_r(i,j)=corr( CityTCSplit{i,j},TempCityRTs);
    end
    
 
    
 %   MountainTC{i,1}=8.85725473293589e-16 + -0.2633*zscore(ind_var_mount(TempMountainItems,1)) + -1.3231*zscore(ind_var_mount(TempMountainItems,2)) + 1.7981*zscore(ind_var_mount(TempMountainItems,3)) + 1.4034*zscore(ind_var_mount(TempMountainItems,4)) + -1.2721*zscore(ind_var_mount(TempMountainItems,5)) + -0.2697*zscore(ind_var_mount(TempMountainItems,6));
 MountainTC=[]; 
 %   Mountain_r(i,1)=nancorr( MountainTC{i,1},TempMountainCEs);
 Mountain_r=[];
 

end
TC_summary.mean_r=mean(City_r,1);
TC_summary.mean_r2z=mean(r2z(City_r),1);
[~,TC_summary.p,TC_summary.CI,TC_summary.tSTATs]=ttest(r2z(City_r));
