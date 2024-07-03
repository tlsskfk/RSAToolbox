function [newParcel,fitMat] = Parcel2ParcelICA(ParcelData,ParcelMask,RefParcelMask,RefStrength)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
RefParcelNums=unique(RefParcelMask(:));
RefParcelNums(RefParcelNums==0,:)=[];
numParcels=length(RefParcelNums);
if istable(ParcelData)
    TCMat=table2array(ParcelData);
else
    TCMat=ParcelData;
end
nanMat=isnan(TCMat);
TPFilter=sum(nanMat,2)==size(nanMat,2);
ParcelFilter=sum(nanMat,1)>size(nanMat,1)*0.5;
TCMat(TPFilter==1,:)=[];
TCMat(:,ParcelFilter==1)=[];

icaOut = rica(TCMat,numParcels,'Standardize',1,'IterationLimit',10000);
icaWeights=nan(length(ParcelFilter),numParcels);
icaWeights(ParcelFilter==0,:)=icaOut.TransformWeights;
[~,newParcelNums]=max(icaWeights,[],2);
SmoothMasks=repmat(RefParcelMask,[1,1,1,numParcels])*0;
for i = 1:numParcels
   tempMask=single(RefParcelMask==RefParcelNums(i,1));
   tempMask(tempMask==1)=100/nansum(tempMask(:));
   SmoothMasks(:,:,:,i)=bvsmoother(tempMask,3);
end
SmoothMasks=reshape(SmoothMasks,[size(SmoothMasks,1)*size(SmoothMasks,2)*size(SmoothMasks,3),size(SmoothMasks,4)]);
newParcel=ParcelMask;
parcelNums=unique(ParcelMask(:));
parcelNums(parcelNums==0,:)=[];
for i = parcelNums'
    newParcel(ParcelMask==i)=newParcelNums(i,1);
end
newParcelVert=newParcel(:);
fitMat=nan(numParcels,numParcels);
for i = 1:numParcels
    fitMat(i,:)=nansum(SmoothMasks(newParcelVert==i,:),1);
end
fitMat_A=fitMat./(repmat(sum(fitMat,2),[1,length(fitMat)]));
fitMat_B=fitMat./(repmat(sum(fitMat,1),[length(fitMat),1]));
fitMatprop=(fitMat_A+fitMat_B)/2;
fitMatVert=[fitMatprop(:),reshape(repmat([1:numParcels],[numParcels,1]),numParcels^2,1),reshape(repmat([1:numParcels]',[1,numParcels]),numParcels^2,1)];
[~,sortRank]=sort(fitMatVert(:,1),1,'descend');
sortedFitMat=fitMatVert(sortRank,:);
newAssignment=nan(numParcels,1);

for i = 1:numParcels
    icaNum=sortedFitMat(1,3);
    parcelNum=sortedFitMat(1,2);
    fitVal=sortedFitMat(1,1);
    newAssignment(icaNum,1)=parcelNum;
    newAssignment(icaNum,2)=fitVal;
    filterMat=(single(sortedFitMat(:,2)==parcelNum)+single(sortedFitMat(:,3)==icaNum))~=0;
    sortedFitMat(filterMat,:)=[];
end
a=1;


end


