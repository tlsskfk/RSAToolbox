function [ VTC_byTrial ] = VTC_byZone(FileName,TRs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(FileName);
VTC_byTrial{1,1}=CPT_analyze_zone_func2(response,data);

temp_median=nanmedian(VTC_byTrial{1,1});
temp_mean=nanmean(VTC_byTrial{1,1});
VTC_byTrial{1,2}=single(VTC_byTrial{1,1}<=temp_mean);
VTC_byTrial{1,3}=single(VTC_byTrial{1,1}<=temp_median);
VTC_byTrial{1,4}=imresize(VTC_byTrial{1,1},[TRs,1]);
temp_mean2=nanmean(VTC_byTrial{1,4});
temp_median2=nanmedian(VTC_byTrial{1,4});
VTC_byTrial{1,5}=single(VTC_byTrial{1,4}<=temp_mean2);
VTC_byTrial{1,6}=single(VTC_byTrial{1,4}<=temp_median2);

end

