function [ sRSMs,xRSMs,Diags,xRSM_Means] = RSMParcer( RSMs,numSplits,sym )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    sym=1;
end


numRSMs=size(RSMs,3);
if sym==0
    nullMat=zeros(size(RSMs,1)*2,size(RSMs,1)*2,numRSMs);
    nullMat(1:size(RSMs,1),size(RSMs,1)+1:end,:)=RSMs;
    RSMs=nullMat;
    nullMat=[];
end


sizeRSM=size(RSMs,1);
numDiags=[0,1,3,6,10,15,21];
numDiags=numDiags(1,numSplits);
numxRSMs=numDiags*2;
numCond=sizeRSM/numSplits;
triLength=(numCond^2-numCond)/2;
sRSMs=zeros(triLength,numSplits,numRSMs,'single');
xRSMs=zeros(triLength,numxRSMs,numRSMs,'single');
Diags=zeros(numCond,numDiags,numRSMs,'single');
xRSM_Means=zeros(numCond,numDiags,numRSMs,'single');

for n=1:numRSMs
    DiagCount=1;
    xRSMCount=1;    
    RSM=RSMs(:,:,n);
    for i = 1:numSplits
        iStart=((i-1)*numCond)+1;
        iEnd=iStart+numCond-1;
        sRSMs(:,i,n)=mat2uppertriuvectormat(RSM(iStart:iEnd,iStart:iEnd));
        if i == numSplits
            continue
        end
        for j = i+1:numSplits
            jStart=((j-1)*numCond)+1;
            jEnd=jStart+numCond-1;
            tempRSM=RSM(iStart:iEnd,jStart:jEnd);
            Diags(:,DiagCount,n)=atanh(tempRSM(eye(numCond)==1));
            tempRSM(eye(numCond)==1)=nan;
            xRSM_Means(:,DiagCount,n)=nanmean(atanh(tempRSM),2);
            DiagCount=DiagCount+1;
            xRSMs(:,xRSMCount,n)=mat2uppertriuvectormat(tempRSM);
            xRSMCount=xRSMCount+1;
            tempRSM=permute(tempRSM,[2,1]);
            xRSMs(:,xRSMCount,n)=mat2uppertriuvectormat(tempRSM);
            xRSMCount=xRSMCount+1;            
        end     
    end
end
if sym==0
    sRSMs=xRSMs;
end
end

