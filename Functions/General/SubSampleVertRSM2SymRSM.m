function [ SymRSMs,vertRSMs ] = SubSampleVertRSM2SymRSM( vertRSM,SubSampleVec,nanDiag )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin==2
    nanDiag=0;
end
RSMsize=length(SubSampleVec);
numRSMs= size(vertRSM,2);
SymRSMs=[];
if isempty(RSMsize)
    N=length(vertRSM(:,1));
    RSMsize=roots([1,-1,-(2*N)]);
    RSMsize(RSMsize<0,:)=[];
end

for i = 1:numRSMs
    SymRSM=ones(int16(RSMsize));
    SymRSM(triu(SymRSM,1)==1)=vertRSM(:,i);
    SymRSM=triRSMs2SymRSMs(SymRSM);
    if nanDiag==1
        SymRSM(eye(length(SymRSM))==1)=nan;
    end
    SymRSM=SymRSM(SubSampleVec==1,SubSampleVec==1);
    SymRSMs=cat(3,SymRSMs,SymRSM);
end
vertRSMs=mat2uppertriuvectormat(SymRSMs);
end

function SymRSMs=triRSMs2SymRSMs(RSMs)
SymRSMs=RSMs;

for i = 1:size(RSMs,3)
    RSM=RSMs(:,:,i);
    SymRSMs(:,:,i)=triu(RSM,1)+triu(RSM)';
end
    
end