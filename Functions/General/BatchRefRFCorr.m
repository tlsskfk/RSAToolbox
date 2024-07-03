function [AllRF] = BatchRefRFCorr(x,y,batchSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
numRF=size(x,2);
AllRF=nan(numRF,1);
batches=1:batchSize:numRF;
batchend=[batches(1,2:end)-1,numRF];
batches=[batches;batchend];
for i = 1:size(batches,2)
    ini=batches(1,i);
    fin=batches(2,i);
    RCAMat=corr(x(:,ini:fin),y(:,ini:fin));
    AllRF(ini:fin,1)=single(RCAMat(eye(length(RCAMat))==1)');
end

end