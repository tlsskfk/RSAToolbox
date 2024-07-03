function [ RSMs ] = ParcellationRSMs_parfor( ActivationMat,ParcellationVector,DistType )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    DistType='corrcoef';
end

ParcellationNums=unique(ParcellationVector);
ParcellationNums(ParcellationNums==0)=[];
numParcels=max(ParcellationNums);
numCond=size(ActivationMat,2);
RSMs=nan(numCond,numCond,numParcels,'single');

if strcmpi(DistType,'corrcoef')
    parfor i = 1:length(ParcellationNums)
        parcelNum=ParcellationNums(i,1);
        RSMs(:,:,i)=corrcoef(ActivationMat(ParcellationVector==parcelNum,:));
    end
else
    ActivationMat=ActivationMat';
    parfor i = 1:length(ParcellationNums)
        parcelNum=ParcellationNums(i,1);
        RSMs(:,:,i)=squareform(pdist(ActivationMat(:,ParcellationVector==parcelNum),DistType));
    end    
end