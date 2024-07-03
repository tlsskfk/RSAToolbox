function [ EventDensity,preEvent,postEvent ] = TargetConcentrationTC( TargetTrials,gauss_mu,gauss_sigma,z_ceiling )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
TargetTrials=TargetTrials(:);
numTarg=sum(TargetTrials(:));
if numTarg > 0
    
gaussWindow = guass_sigma * 4;
LeftGauss = normpdf([gaussWindow*-1:0],gauss_mu,gauss_sigma);
RightGauss = normpdf([0:gaussWindow],gauss_mu,gauss_sigma);
FullGauss = normpdf([gaussWindow*-1:gaussWindow],gauss_mu,gauss_sigma);

else
    EventDensity=zeros(length(TargetTrials),1);
    preEvent=zeros(length(TargetTrials),1);
    postEvent=zeros(length(TargetTrials),1);    
end


c=0;
new=zeros(numNonTarg+1,1);
count=1;
for i = 1: size(TargetTrials,1)
    m=TargetTrials(i,1);
    if m==1
        c=c+1;
    else
        new(count,1)=c;
        c=0;
        count=count+1;
    end;
end;
final=zeros(numNonTarg+1,numNonTarg+1);
for i= 1:numNonTarg+1
    xxx=0:1:(numNonTarg+1-i);
    gauss=normpdf(xxx,gauss_mu,gauss_sigma)*new(i,1);
    final(i,i:end)=gauss;
end;

final2=sum(final);
final2(end)=[];
TargetFreq_RT_TC=zscore(final2);
TargetFreq_RT_TC(TargetFreq_RT_TC>z_ceiling)=z_ceiling;
Timecourse=TargetFreq_RT_TC;
end

