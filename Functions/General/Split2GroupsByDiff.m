function [ Splits,GroupMeans,RealDiff,Error ] = Split2GroupsByDiff( Vals,TargetDiff,subSampleSize,reps,ErrorThresh )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    reps=10000;
    ErrorThresh=0;
    subSampleSize=0.5;
end
if nargin==3
    ErrorThresh=0;
    reps=10000;
end
if nargin==4
    ErrorThresh=0;
end


if subSampleSize<1
    subSampleSize=subSampleSize*100;
end

numVals=length(Vals);
iniRand=rand(numVals,1);
if rand(1) > 0.5
    iGroup = iniRand < prctile(iniRand,subSampleSize);
else
    iGroup = iniRand <= prctile(iniRand,subSampleSize);
end
maxFlip=3;
if subSampleSize >= 50
    if sum(single(iGroup==0)) < 3
        maxFlip=sum(single(iGroup==0));
    end
else
    if sum(single(iGroup==1)) < 3
        maxFlip=sum(single(iGroup==1));
    end
end
iDiff=nanmean(Vals(iGroup==1),1)-nanmean(Vals(iGroup==0),1);
iError=abs(iDiff-TargetDiff);
for rep=1:reps
    A=find(iGroup==0);
    B=find(iGroup==1);
    fGroup=iGroup;
    count=0;
    while count < randi(maxFlip)
        fGroup(A(randi(length(A))),1)=1;
        fGroup(B(randi(length(B))),1)=0;
        A=find(fGroup==0);
        B=find(fGroup==1);
        count=count+1;
    end
    fDiff=nanmean(Vals(fGroup==1),1)-nanmean(Vals(fGroup==0),1);
    fError=abs(fDiff-TargetDiff);
    if fError<iError
       iError=fError; 
       iGroup=fGroup;
       iDiff=fDiff;
    end
    if iError<=ErrorThresh
        break;
    end
end
Splits=iGroup;
RealDiff=iDiff;
Error=iError;
GroupMeans=[nanmean(Vals(iGroup==1),1),nanmean(Vals(iGroup==0),1)];

