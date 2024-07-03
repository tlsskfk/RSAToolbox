function [ SigMat ] = SigSeg2( pMat,pThreshs,signType)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SigMat=pMat*0;
for i = 1:length(pThreshs)
    if strcmpi(signType,'both')
        SigMat=SigMat+single(pMat<=pThreshs(1,i)|pMat>=1-pThreshs(1,i));
    elseif strcmpi(signType,'pos')
        SigMat=SigMat+single(pMat<=pThreshs(1,i));
    elseif strcmpi(signType,'neg')
        SigMat=SigMat+single(pMat>=1-pThreshs(1,i));
    end
end
% 
% if strcmpi(signType,'both')
%     SigMat=SigMat.*single(sign(rMat));
% elseif strcmpi(signType,'pos')
%     SigMat=SigMat.*single(sign(rMat)==1);
% elseif strcmpi(signType,'neg')
%     SigMat=SigMat.*single(sign(rMat)==-1);
% end