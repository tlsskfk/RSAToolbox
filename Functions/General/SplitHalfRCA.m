function [RF,RC,IF,RFs,RFx,RCs,RCx] = SplitHalfRCA(RSMs)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ sRSMs,xRSMs,Diags,xRSM_Means] = RSMParcer( RSMs,2 );
sRSMs=permute(sRSMs,[1,3,2]);
[ ~,RCs,RFs,~] = RSMParcer( corrcoef([sRSMs(:,:,1),sRSMs(:,:,2)]),2);
RFs=atanh(RFs(:,1));
RCs=nanmean(atanh(RCs),2);
xRSMs=permute(xRSMs,[1,3,2]);
[ ~,RCx,RFx,~] = RSMParcer( corrcoef([xRSMs(:,:,1),xRSMs(:,:,2)]),2);
RFx=atanh(RFx(:,1));
RCx=nanmean(atanh(RCx),2);
IF=squeeze(nanmean(Diags-xRSM_Means,1));
% xRSM_Means=squeeze(nanmean(xRSM_Means,1));
% IF=Diags-xRSM_Means;
RF=nanmean([RFs,RFx],2);
RC=nanmean([RCs,RCx],2);
end

