function [meanSubSampleVar,maxSubSampleVar] = Compute_SubSampleVariable(InVar,SubSamples)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if iscell(InVar)
    tempInVar=InVar;
    InVar=[];
    for i = 1:size(tempInVar,1)
        if istable(tempInVar{i,1})
            InVar=[InVar;table2array(tempInVar{i,1})];
        else
            InVar=[InVar;tempInVar{i,1}];
        end
    end
end    
if iscell(SubSamples)
    tempSubSamples=SubSamples;
    SubSamples=[];
    for i = 1:size(tempSubSamples,1)
        if istable(tempSubSamples{i,1})
            SubSamples=[SubSamples;table2array(tempSubSamples{i,1})];
        else
            SubSamples=[SubSamples;tempSubSamples{i,1}];
        end
    end
end 

numVar=size(InVar,2);
numSubSamples=size(SubSamples,2);
meanSubSampleVar=nan(numSubSamples,numVar);
maxSubSampleVar=nan(numSubSamples,numVar);
for i = 1:numSubSamples
    meanSubSampleVar(i,:)=nanmean(InVar(SubSamples(:,i)==1,:),1);
    maxSubSampleVar(i,:)=max(InVar(SubSamples(:,i)==1,:));
end    
end

