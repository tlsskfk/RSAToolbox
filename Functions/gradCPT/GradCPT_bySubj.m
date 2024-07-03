function [Results] = GradCPT_bySubj( FileName,CitySimVals,MountainSimVals,split,VTCs,RT_TCs,CE_TCs,ComputeFullVar )
if nargin < 4
    ComputeFullVar=0;
    VTCs=[];
    split=[];
    RT_TCs=[];
    CE_TCs=[];
end
if nargin < 5
    ComputeFullVar=0;
    VTCs=[];
    RT_TCs=[];
    CE_TCs=[];
end
if nargin < 6
    ComputeFullVar=0;
    RT_TCs=[];
    CE_TCs=[];
end
if nargin < 7
    ComputeFullVar=0;
    CE_TCs=[];
end
if nargin < 8
    ComputeFullVar=0;
end
if isstruct(FileName)
    data=FileName.data;
    response=FileName.response;
    ttt=FileName.ttt;
    if isfield(FileName,'bordertracker')
        bordertracker=FileName.bordertracker;
    end
else
    load(FileName);
end
%[ bins ] = makebins( size(data,1),split(1,2) );
if ~isempty(split)
    data=data(split==1,:);
    response=response(split==1,:);
    ttt=ttt(split==1,:);
    if exist('bordertracker','var')
        if length(bordertracker)==length(split)
            Reward=bordertracker(:,2)==255; 
        else
            Reward=zeros(length(split),1);
        end
    else
        Reward=zeros(length(split),1);
    end
    Reward=Reward(split==1,:);
    if ~isempty(VTCs)
        VTCs=VTCs(split==1,:);
    end
else
    if exist('bordertracker','var')
        if length(bordertracker)==length(ttt)
            Reward=bordertracker(:,2)==255; 
        else
            Reward=zeros(length(ttt),1);
        end
    else
        Reward=zeros(length(ttt),1);
    end
end
Results.RewardSum=sum(single(Reward(:)));
Results.RewardTC=Reward;  
StimItems=data(:,4);
Results.CityItems = StimItems(data(:,3)==2,1);
Results.MountainItems = StimItems(data(:,3)==1,1);
Results.RT_byTrial = response(:,5)*1000;
Results.RT_byTrial(Results.RT_byTrial==0) = nan;

Results.MountainTrials=data(:,3)==1;
Results.CityTrials=data(:,3)==2;

Results.CityRTs=Results.RT_byTrial(response(:,1)==2);

if isempty(RT_TCs)==0
    tempReg=regstats(Results.CityRTs,RT_TCs,'linear',{'r','beta'});
    Results.CityRTs=tempReg.r+tempReg.beta(1,1);
    response(response(:,1)==2,5)=Results.CityRTs/1000;
    response(isnan(response(:,5)),5)=0;
end

Results.CityRTs_Z = nan_zscore(Results.CityRTs);
Results.MountainRTs=Results.RT_byTrial(response(:,1)==1);
Results.MountainRTs_Z=nan_zscore(Results.MountainRTs);
Results.CityOE=single(isnan(Results.CityRTs));
Results.MountainCE=single(isnan(Results.MountainRTs)==0);

if isempty(CE_TCs)==0
    tempReg=regstats(Results.MountainCE,CE_TCs,'linear',{'r','beta'});
    Results.MountainCE=tempReg.r;
end

Results.CityRT_byItem=zeros(10,1);
Results.MountainRT_byItem=zeros(10,1);
Results.CityRT_byItem_Z=zeros(10,1);
Results.MountainRT_byItem_Z=zeros(10,1);
Results.CityOE_byItem=zeros(10,1);
Results.MountainCE_byItem=zeros(10,1);
Results.CityN_byItem=zeros(10,1);
Results.MountainN_byItem=zeros(10,1);

for i =1:10
    Results.CityRT_byItem(i,1)=nanmean(Results.CityRTs(Results.CityItems==i));
    Results.MountainRT_byItem(i,1)=nanmean(Results.MountainRTs(Results.MountainItems==i));
    Results.CityRT_byItem_Z(i,1)=nanmean(Results.CityRTs_Z(Results.CityItems==i));
    Results.MountainRT_byItem_Z(i,1)=nanmean(Results.MountainRTs_Z(Results.MountainItems==i));
    Results.CityOE_byItem(i,1)=nansum(Results.CityOE(Results.CityItems==i))/size(Results.CityOE(Results.CityItems==i),1);
    %Results.MountainCE_byItem(i,1)=nansum(Results.MountainCE(Results.MountainItems==i))/size(Results.MountainCE(Results.MountainItems==i),1);
    Results.MountainCE_byItem(i,1)=nansum(Results.MountainCE(Results.MountainItems==i));
    Results.CityN_byItem(i,1)=size(Results.CityOE(Results.CityItems==i),1);
    Results.MountainN_byItem(i,1)=size(Results.MountainCE(Results.MountainItems==i),1);
end
if isempty(VTCs)
    [Results.VTC_byTrial,~,medVar]=CPT_analyze_zone_func2(response,data);
    Results.VTC_byTrial=imresize(Results.VTC_byTrial,[size(data,1),1]);
else
    [~,~,medVar]=CPT_analyze_zone_func2(response,data);
    Results.VTC_byTrial=VTCs;
end
Results.CityVTC=Results.VTC_byTrial(response(:,1)==2);
Results.MountainVTC=Results.VTC_byTrial(response(:,1)==1);
% Results.VTC_byTrial=[];
% Results.CityVTC=[];
% Results.MountainVTC=[];
[Results.Zone,Results.ZonebyItem,Results.ZonePrime,Results.ZonePrimebyItem,Results.VTCderiv] = ZoneByItem( Results.VTC_byTrial,Results.CityTrials,StimItems );
VTC_IN=median(Results.VTC_byTrial(Results.Zone==1));
VTC_OUT=median(Results.VTC_byTrial(Results.Zone==0));
VTCzonediff=VTC_IN-VTC_OUT;
VTC_itemIN=median(Results.VTC_byTrial(Results.ZonebyItem==1));
VTC_itemOUT=median(Results.VTC_byTrial(Results.ZonebyItem==0));
VTCitemzonediff=VTC_itemIN-VTC_itemOUT;

VTCderiv_WAX=median(Results.VTCderiv(Results.ZonePrime==1));
VTCderiv_WAN=median(Results.VTCderiv(Results.ZonePrime==0));
VTCderivPhaseDiff=VTCderiv_WAX-VTCderiv_WAN;
VTCderiv_itemWAX=median(Results.VTCderiv(Results.ZonePrimebyItem==1));
VTCderiv_itemWAN=median(Results.VTCderiv(Results.ZonePrimebyItem==0));
VTCderivItemPhaseDiff=VTCderiv_itemWAX-VTCderiv_itemWAN;


Var_byItem_city=[Results.CityRT_byItem];
Var_byItem_mountain=[Results.MountainCE_byItem];
Var_byTrial_city = [Results.CityRTs];
Var_byTrial_mountain = [Results.MountainCE];

[Results.byItem_city_r,Results.byItem_city_p]=corr(Var_byItem_city,CitySimVals,'rows','pairwise');
%[Results.byItem_mountain_r,byItem_mountain_p]=nancorr(Var_byItem_mountain,MountainSimVals);
Results.byItem_mountain_r=[];
Results.byItem_mountain_p=[];
[Results.byTrial_city_r,Results.byTrial_city_p]=corr(Var_byTrial_city,CitySimVals(Results.CityItems,:),'rows','pairwise');
if sum(Var_byTrial_mountain(:))<2
    Results.byTrial_mountain_r=nan(1,size(MountainSimVals,2));
    Results.byTrial_mountain_p=nan(1,size(MountainSimVals,2));
else
    [Results.byTrial_mountain_r,Results.byTrial_mountain_p]=corr(Var_byTrial_mountain,MountainSimVals(Results.MountainItems,:),'rows','pairwise');
end



[Output1]=CPT_analyze_func2(response,ttt);
if ComputeFullVar==1
    [ outputIN ] = Compute_Standard_gradCPT_Vars(response(Results.Zone==1,:));
    [ outputOUT ] = Compute_Standard_gradCPT_Vars(response(Results.Zone==0,:));
    [ output_itemIN ] = Compute_Standard_gradCPT_Vars(response(Results.ZonebyItem==1,:));
    [ output_itemOUT ] = Compute_Standard_gradCPT_Vars(response(Results.ZonebyItem==0,:));
    [ outputWAX ] = Compute_Standard_gradCPT_Vars(response(Results.ZonePrime==1,:));
    [ outputWAN ] = Compute_Standard_gradCPT_Vars(response(Results.ZonePrime==0,:));
    [ output_itemWAX ] = Compute_Standard_gradCPT_Vars(response(Results.ZonePrimebyItem==1,:));
    [ output_itemWAN ] = Compute_Standard_gradCPT_Vars(response(Results.ZonePrimebyItem==0,:));
    Results.Output_Full=[Output1,medVar,outputIN,VTC_IN,outputOUT,VTC_OUT,outputIN-outputOUT,VTCzonediff,output_itemIN,VTC_itemIN,output_itemOUT,VTC_itemOUT,output_itemIN-output_itemOUT,VTCitemzonediff,median(Results.VTCderiv),outputWAX,VTCderiv_WAX,outputWAN,VTCderiv_WAN,outputWAX-outputWAN,VTCderivPhaseDiff,output_itemWAX,VTCderiv_itemWAX,output_itemWAN,VTCderiv_itemWAN,output_itemWAX-output_itemWAN,VTCderivItemPhaseDiff];
else
    Results.Output_Full=[Output1,medVar,VTC_IN,VTC_OUT,VTCzonediff,VTC_itemIN,VTC_itemOUT,VTCitemzonediff,median(Results.VTCderiv),VTCderiv_WAX,VTCderiv_WAN,VTCderivPhaseDiff,VTCderiv_itemWAX,VTCderiv_itemWAN,VTCderivItemPhaseDiff];    
end
Results.Output=Results.Output_Full(1,[9,14,19,24,29,34,39,44]);
[Results.ZO_count] = ZOs(response,38,0.01);
%Output=[Pre_RTs Post_RTs commission_rate omission_rate meanRT STD_RT error_rate dprime criterion CV PES CE_Slope OE_Slope RT_Slope STD_Slope Err_Slope dprime_Slope criterion_Slope CV_Slope];

end

