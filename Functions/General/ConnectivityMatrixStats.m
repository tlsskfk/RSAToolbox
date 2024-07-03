function [ResultsByNode,ResultsByCM] = ConnectivityMatrixStats(CM,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ResultsByNode=[];
if size(CM,1)~=size(CM,2)
    CM=vertRSM2SymRSM(CM);
end
signMat=double(sign(CM));
[pmat] = VariableSetter('pmat',[],varargin);
if ~isempty(pmat) && size(pmat,1)~=size(pmat,2)
    pmat=vertRSM2SymRSM(pmat);
end 
CM(eye(length(CM))==1)=nan;
MaxDegreeByNode=size(CM,1)-1;
MaxDegreeByCM=size(CM,1)^2 - size(CM,1);
if ~isempty(pmat)
    [ThresholdType] = VariableSetter('ThresholdType','pval',varargin); %stat or rank, pval
    [Threshold] = VariableSetter('Threshold',0.05,varargin);
    binCM=pmat<=Threshold; 
    signedCM=double(binCM).*signMat;
    posCM=signedCM==1;
    negCM=signedCM==-1;  
else
    %ThresholdType
    [ThresholdType] = VariableSetter('ThresholdType','rank',varargin); %stat, or rank, pval
    %Threshold 
    if strcmpi(ThresholdType,'rank')
        [Threshold] = VariableSetter('Threshold',0.10,varargin);
        Threshold=Threshold*100;
        ThreshVals=prctile(CM(:),[Threshold,100-Threshold]);
        negCM=CM<=ThreshVals(1,1);
        posCM=CM>=ThreshVals(1,2);
        binCM=(single(negCM)+single(posCM))==1;
        signedCM=double(binCM).*signMat;
    else
        [Threshold] = VariableSetter('Threshold',[-0.10,0.10],varargin);
        negCM=CM<=Threshold(1,1);
        posCM=CM>=Threshold(1,2);
        binCM=(single(negCM)+single(posCM))==1;
        signedCM=double(binCM).*signMat;        
    end        
end



ResultsByNode.PosDeg=nansum(single(posCM),2)/MaxDegreeByNode;
ResultsByNode.NegDeg=nansum(single(negCM),2)/MaxDegreeByNode;
ResultsByNode.TotDeg=nansum(single(binCM),2)/MaxDegreeByNode;
ResultsByNode.NetDeg=nansum(single(signedCM),2)/MaxDegreeByNode;
ResultsByNode.DepDeg=(nansum(single(binCM),2)-nansum(single(signedCM),2))/MaxDegreeByNode;
ResultsByNode.PosMean=nansum(single(posCM).*CM,2)./ResultsByNode.PosDeg;
ResultsByNode.PosMean(isnan(ResultsByNode.PosMean))=0;
ResultsByNode.NegMean=nansum(single(negCM).*CM,2)./ResultsByNode.NegDeg;
ResultsByNode.NegMean(isnan(ResultsByNode.NegMean))=0;
ResultsByNode.NetMean=nanmean(CM,2);
ResultsByNode.NetMean(isnan(ResultsByNode.NetMean))=0;
ResultsByNode=struct2table(ResultsByNode);

ResultsByCM.PosDeg=nansum(single(posCM(:)))/MaxDegreeByCM;
ResultsByCM.NegDeg=nansum(single(negCM(:)))/MaxDegreeByCM;
ResultsByCM.TotDeg=nansum(single(binCM(:)))/MaxDegreeByCM;
ResultsByCM.NetDeg=nansum(single(signedCM(:)))/MaxDegreeByCM;
ResultsByCM.DepDeg=(nansum(single(binCM(:)))-nansum(single(signedCM(:))))/MaxDegreeByCM;
ResultsByCM.PosMean=nansum(single(posCM(:)).*CM(:))./ResultsByCM.PosDeg;
ResultsByCM.PosMean(isnan(ResultsByCM.PosMean))=0;
ResultsByCM.NegMean=nansum(single(negCM(:)).*CM(:))./ResultsByCM.NegDeg;
ResultsByCM.NegMean(isnan(ResultsByCM.NegMean))=0;
ResultsByCM.NetMean=nanmean(CM(:));
ResultsByCM.NetMean(isnan(ResultsByCM.NetMean))=0;
ResultsByCM=struct2table(ResultsByCM);



end

