function [AllMeans,AllZs] = SubSampleCondRF(RSMs,SubSize,Reps)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if iscell(RSMs)
    RSMs=cell2nDMAT(RSMs);
end
numConds=size(RSMs,1);
numROIs=size(RSMs,3);
numSS=size(RSMs,4);
AllMeans=nan(Reps,numROIs);
AllZs=nan(Reps,numROIs);
parfor rep = 1:Reps
    ind=randperm(numConds,SubSize);
    [ ~,RF] = FastRCA(RSMs(ind,ind,:,:),'vertIn',0);
    [~,AllZs(rep,:),AllMeans(rep,:)]=getTval(RF',1)    
end
AllMeans=num2cell(AllMeans,1);
AllZs=num2cell(AllZs,1);
end

