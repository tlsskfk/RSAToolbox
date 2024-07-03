function [ RSMs ] = ParcellationRSMs( ActivationMat,ParcellationVector,DistType,numParcels )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ParcellationNums=unique(ParcellationVector);
ParcellationNums(ParcellationNums==0)=[];
if nargin==2
    DistType='corrcoef';
    numParcels=max(ParcellationNums);
end
if nargin==3
    numParcels=max(ParcellationNums);
end
numCond=size(ActivationMat,2);
RSMs=nan(numCond,numCond,numParcels,'single');
if ~isempty(ParcellationNums)
    if strcmpi(DistType,'corrcoef')
        for i = ParcellationNums'
            RSMs(:,:,i)=corrcoef(ActivationMat(ParcellationVector==i,:));
        end
    elseif strcmpi(DistType,'meanSim')
        for i = ParcellationNums'
            RSMs(:,:,i)=vector2diffdistmat(nanmean(ActivationMat(ParcellationVector==i,:),1));
        end        
    else    
        ActivationMat=ActivationMat';
        for i = ParcellationNums'
            RSMs(:,:,i)=squareform(pdist(ActivationMat(:,ParcellationVector==i),DistType));
        end    
    end
end
end