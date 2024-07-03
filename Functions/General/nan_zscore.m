function [ zmat ] = nan_zscore( mat,type )
%Written by David Rothlein
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
  type = '';
end
i=size(mat,1);

if strcmp(type,'pooled')==1
    pooled_mean=nanmean(mat(:));
    pooled_std=nanstd(mat(:));
    zmat=(mat-pooled_mean)/pooled_std;
else
    nanmeans=repmat(nanmean(mat,1),i,1);
    nanstds=repmat(nanstd(mat),i,1);    
    zmat=(mat-nanmeans)./nanstds;
end
end

