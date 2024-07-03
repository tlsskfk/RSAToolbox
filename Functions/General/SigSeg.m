function [ SigMat,SigLabelMat ] = SigSeg( pMat,pThreshs,signType,ThreshLabels)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    pThreshs=[0.05,0.01,0.005,0.001];
    signType='pos';
    for i = 1:length(pThreshs)
        ThreshLabels{1,i}=repmat('*',[1,i]);
    end
end
if nargin==2
    signType='pos';
    for i = 1:length(pThreshs)
        ThreshLabels{1,i}=repmat('*',[1,i]);
    end
end
if nargin==3
    for i = 1:length(pThreshs)
        ThreshLabels{1,i}=repmat('*',[1,i]);
    end
end

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

SigLabelMat=cell(size(SigMat));
for i =1:size(SigMat,1)
    for j = 1:size(SigMat,2)
        if SigMat(i,j)>0
            SigLabelMat{i,j}=ThreshLabels{1,SigMat(i,j)};
        else
            SigLabelMat{i,j}=[];
        end
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