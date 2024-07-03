function [AlignedRSMs,AlignInds,OverlapMask] = Searchlight_RSMAlign(RSMs,RSM_masks)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
OverlapMask=single(sum(RSM_masks,4)>0);
emptyRSMInd=find(~cellfun(@isempty,RSMs));
AlignedRSMs=nan(size(RSMs{emptyRSMInd(1,1),1},1),sum(OverlapMask(:)),size(RSMs,1),'single');
AlignInds=zeros(sum(OverlapMask(:)),size(RSMs,1));
for f=1:size(RSMs,1)
    if ~isempty(RSMs{f,1})
        tempMask=RSM_masks(:,:,:,f);
        AlignInd=tempMask(OverlapMask==1);
        AlignedRSMs(:,AlignInd==1,f)=RSMs{f,1};
        AlignInds(:,f)=AlignInd;
    end
end
end