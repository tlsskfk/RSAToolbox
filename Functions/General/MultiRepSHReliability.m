function [RelVal,AllRs] = MultiRepSHReliability(InData,UseCorrection)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    UseCorrection=0;
end
numItems=size(InData,2);
SplitVec=ones(1,numItems);
HalfInd=ceil(numItems/2);
SplitVec(1,1:HalfInd)=2;
numReps=1000;
AllRs=zeros(numReps,1);
for rep=1:numReps
   RandVec=SplitVec(1,randperm(length(SplitVec)));
   a=nanmean(InData(:,RandVec==1),2);
   b=nanmean(InData(:,RandVec==2),2);
   r=corr(a,b);
   if UseCorrection == 1
        AllRs(rep,1)=(2*r)/(1+r);
   else
        AllRs(rep,1)=r;
   end
end
RelVal=mean(AllRs,1);

end

