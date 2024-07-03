function [AllboldTCs,commonMask] = GLM_boldTCPrep(boldTCs,boldMasks,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ParcellationMask] = VariableSetter('ParcellationMask',[],varargin);
[ResampleSizes] = VariableSetter('ResampleSizes',[],varargin);
[NormWithinRun] = VariableSetter('NormWithinRun','',varargin);
[NormAcrossRun] = VariableSetter('NormAcrossRun','',varargin);
% Remove N volumes at the beginning of run
[RemoveStartVols] = VariableSetter('RemoveStartVols',[],varargin);

if ~iscell(boldTCs)
    boldTCs={boldTCs};
    boldMasks={boldMasks};
end

numRuns=size(boldTCs,1);
for run=1:numRuns
    if max((boldMasks{run,1}(:)))==1
        boldMasks{run,1}(boldMasks{run,1}==1)=[1:sum(boldMasks{run,1}(:))];
    end
end

commonMask=single(sum(single(cell2nDMAT(boldMasks)~=0),4)==numRuns);
AllboldTCs=[];
if ~isempty(ParcellationMask)
    commonMask=commonMask.*single(ParcellationMask~=0);
end

for run=1:numRuns
    useInd=boldMasks{run,1}(find(commonMask));
    useInd=useInd(useInd~=0);
    tempboldTC=boldTCs{run,1}(:,useInd);
    if ~isempty(RemoveStartVols)
        tempboldTC(1:RemoveStartVols,:)=nan;
    end
    if strcmpi(NormWithinRun,'zscore')
        tempboldTC=nan_zscore(tempboldTC,1);
    end
    if iscell(ResampleSizes)
        tempboldTC=imresize(tempboldTC,[ResampleSizes{run,1},size(tempboldTC,2)]);
    end
    AllboldTCs=[AllboldTCs;tempboldTC];
end 

if strcmpi(NormAcrossRun,'zscore')
    AllboldTCs=nan_zscore(AllboldTCs,1);
end
commonMask(commonMask==1)=[1:sum(commonMask(:))];    
end

