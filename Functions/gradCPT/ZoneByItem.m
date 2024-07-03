function [Zone,ZonebyItem,ZonePrime,ZonePrimebyItem,VTCderiv] = ZoneByItem( VTC,StimType,StimNums,setDiff )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    setDiff=[];
end

if length(VTC) < length(StimType)
    toRemove=length(StimType)-length(VTC)-1;
    StimType(end-toRemove:end,:)=[];
    StimNums(end-toRemove:end,:)=[];
end

Zone = VTC<nanmedian(VTC);
VTCshift=[VTC(2:end,1);median(VTC)];
VTCderiv=VTCshift-VTC; %(VTC(t+1) - VTC(t))/1t or the change in VTC over the change in trials for each trial
ZonePrime=VTCderiv<nanmedian(VTCderiv); %1 = Waxing (Decreasing variability) %0 = Waning (Increasing Vraiability)
numstim=max(StimNums);
ZonebyItem=VTC*0;
ZonePrimebyItem=VTC*0;
for stimnum=1:numstim
    TempZone=VTC(StimType==2 & StimNums==stimnum,:);
    TempZonePrime=VTCderiv(StimType==2 & StimNums==stimnum,:);
    if isempty(setDiff)
        ZonebyItem=ZonebyItem+single(VTC<nanmedian(TempZone) & StimType==2 & StimNums==stimnum);
        ZonePrimebyItem=ZonePrimebyItem+single(VTCderiv<nanmedian(TempZonePrime) & StimType==2 & StimNums==stimnum);
    else
        [ GroupAssigns,~,~ ] = FindGroupDiff(TempZone,setDiff,1000);
        GroupAssigns=(((GroupAssigns*2-1)*-1)+1)/2; % flip assivalue assignment so 1 = In Zone and 0 = out Zone
        ZonebyItem(StimType==2 & StimNums==stimnum,:)=GroupAssigns;
%         [ GroupAssigns,~,~ ] = FindGroupDiff(TempZonePrime,setDiff,1000);
%         GroupAssigns=(((GroupAssigns*2-1)*-1)+1)/2; % flip assivalue assignment so 1 = In Zone and 0 = out Zone        
%         ZonePrimebyItem(StimType==2 & StimNums==stimnum,:)=GroupAssigns;
    end

end    


