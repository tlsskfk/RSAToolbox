function [CorrVals,CorrValsVert,meanVals,zVals] = ComputeFC(TimeSeries,varargin)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

getZVal = VariableSetter( 'getZVal',0,varargin);
VolumeSelect = VariableSetter( 'VolumeSelect',[],varargin);
Lags = VariableSetter( 'Lags',0,varargin);
ms2sec = VariableSetter( 'ms2sec',0,varargin);
SkipFirstN = VariableSetter( 'SkipFirstN',0,varargin);

TimeSeries(isinf(TimeSeries))=nan;
numROIs=size(TimeSeries,2);
numLags=length(Lags);
CorrVals=nan(numROIs,numROIs,numLags,'single');
CorrValsVert=nan((numROIs^2-numROIs)/2,numLags,'single');
meanVals=nan(numROIs,numLags,'single');
zVals=nan(numROIs,numLags,'single');
if ~isempty(VolumeSelect)
    if ms2sec==100
        VolumeSelect=imresize(VolumeSelect,[round(length(VolumeSelect)/10),1],'nearest');
    end
    if size(TimeSeries,1) ~= size(VolumeSelect,1)
        TimeSeries = imresize(TimeSeries,[size(VolumeSelect,1),numROIs]);
    end
    if SkipFirstN~=0
        VolumeSelect(1:SkipFirstN,:)=[];
        TimeSeries(1:SkipFirstN,:)=[];
    end
    VolumeSelectInd=find(VolumeSelect==1);
    for lagNum = 1:numLags
        tempVolumeSelectInd=VolumeSelectInd+Lags(1,lagNum);
        tempVolumeSelectInd(tempVolumeSelectInd<1 | tempVolumeSelectInd>=length(VolumeSelect),:)=[];
        CorrInput=TimeSeries(tempVolumeSelectInd,:);
        CorrInput(isinf(CorrInput))=nan;
        if sum(isnan(CorrInput(:)))>0
            tempCorr=corrcoef(CorrInput,'rows','pairwise');
            meanVals(:,lagNum)=nanmean(CorrInput,1)';
        else
            tempCorr=corrcoef(CorrInput);
            meanVals(:,lagNum)=mean(CorrInput,1)';
        end
        if getZVal==1
            [~,Z]=getTval(CorrInput,1);            
            zVals(:,lagNum)=Z';
        end            
        CorrVals(:,:,lagNum)=tempCorr;
        CorrValsVert(:,lagNum)=mat2uppertriuvectormat(tempCorr);        
    end
else
    if SkipFirstN~=0
        TimeSeries(1:SkipFirstN,:)=[];
    end
    if sum(isnan(TimeSeries(:)))>0
        CorrVals=corrcoef(TimeSeries,'rows','pairwise');
        meanVals(:,1)=nanmean(TimeSeries,1)';
    else
        CorrVals=corrcoef(TimeSeries);
        meanVals(:,1)=mean(TimeSeries,1)';
    end
    CorrVals(eye(length(CorrVals))==1)=nan;
CorrValsVert=mat2uppertriuvectormat(CorrVals);  
end    
