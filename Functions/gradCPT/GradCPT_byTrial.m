function [ CityItems,MountainItems,RT_byTrial,CityRTs,MountainRTs,CityOE,MountainCE,VTC_byTrial,CityVTC,MountainVTC] = GradCPT_byTrial( FileName )

load(FileName);
CityItems = data(data(:,3)==2,4);
MountainItems = data(data(:,3)==1,4);
RT_byTrial = response(:,5)*1000;
RT_byTrial(RT_byTrial==0) = nan;
CityRTs=RT_byTrial(response(:,1)==2);
MountainRTs=RT_byTrial(response(:,1)==1);
CityOE=single(isnan(CityRTs));
MountainCE=single(isnan(MountainRTs)==0);

VTC_byTrial=CPT_analyze_zone_func2(response,data);
CityVTC=VTC_byTrial(response(:,1)==2);
MountainVTC=VTC_byTrial(response(:,1)==1);
end

