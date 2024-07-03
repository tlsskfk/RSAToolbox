function [NN,MeanIJ,MinIJ,MaxIJ,absDiffMean,absDiffMeanInv]=Compute1D_RSMs(VarMat,ssNames)
%Written by David Rothlein
    % NN = similar scores are similar
    % MeanIJ = scores are high for the same reason (high scores are similar)
    % since mean IJ is symmetric, neg RSM corr can be interpretted as opposite.
    % MinIJ = scores are high for the same reason but the lowest score sets
    % the similarity (only as good as the weakest link!)
   
    % absDiffMean = scores are high for the same reason and only similar
    % scores are similar. (i.e. scores must be similar AND both scores must
    % be high)
    % absDiffMeanInv = scores are low for the same reason and only similar
    % scores are similar. (i.e. scores must be similar AND both scores must
    % be low)    
    [~,~,~,Nodes]=labels2uppertriuvectorlabels( ssNames );
    VarMat=scaleVals(VarMat,1,0,1);
    numVars=size(VarMat,2);
    NN=nan(size(Nodes,1),numVars,'single'); 
    MeanIJ=nan(size(Nodes,1),numVars,'single');
    MinIJ=nan(size(Nodes,1),numVars,'single');
    MaxIJ=nan(size(Nodes,1),numVars,'single'); 
    for i = 1:numVars
        tempVar=VarMat(:,i);
        v1=tempVar(Nodes(:,1),1);
        v2=tempVar(Nodes(:,2),1);        
        NN(:,i) = 1 - abs(v1-v2);
        MeanIJ(:,i) = nanmean([v1,v2],2);
        MinIJ(:,i) = min([v1,v2],[],2);
        MaxIJ(:,i) = 1 - max([v1,v2],[],2);
    end
    absDiffMean=NN.*MeanIJ;
    absDiffMeanInv=NN.*(1-MeanIJ);
end