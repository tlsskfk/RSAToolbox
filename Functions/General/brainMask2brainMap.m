function [brainMap] = brainMask2brainMap(brain_vals,brain_mask)
%Written by David Rothlein
%UNTITLED Summary of this function goes here

%   Detailed explanation goes here
numConds=size(brain_vals,2);
brainMap = repmat(brain_mask*0,[1,1,1,numConds]);
for i = 1:numConds
    tempMat=brain_mask;
    tempMat(tempMat==0)=nan;
    tempMat(~isnan(tempMat)) = brain_vals(:,i);
    brainMap(:,:,:,i)=tempMat;
end
end

