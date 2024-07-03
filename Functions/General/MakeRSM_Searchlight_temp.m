function [ RSMs,rsm_mask,numVoxels] = MakeRSM_Searchlight_temp( vals,Mask,SLrad,SLThresh,Shape,MaxBatchSize)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%BatchSize of 5 = fastest
%vals= activation patterns (Voxel,Condition)
numConds=size(vals,2);
RSMsize=(numConds^2-numConds)/2;
%Generate searchlight indicies from mask
[SLinds,numVoxels,MaskCoords] = Mask2Searchlight(Mask,SLrad,SLThresh,Shape);
CleanInd = numVoxels < SLThresh;
SLinds(CleanInd,:)=[];
MaskCoords(CleanInd,:)=[];
numVoxels(CleanInd,:)=[];
numSLs=length(SLinds);
SLsizes=unique(numVoxels);
if useGPU==1
    gpuDevice;
    RSMs=zeros(RSMsize,numSLs,'gpuArray');
else
    RSMs=zeros(RSMsize,numSLs,'single');
end
rsm_mask=coords2mat(MaskCoords,Mask,ones(numSLs,1));
for i = SLsizes'
    selectInd=find(numVoxels==i);
    runSize=length(selectInd);       
    if runSize<=MaxBatchSize
        tempSLinds=SLinds(selectInd,:);
        actMat=zeros(i,runSize*10,'single');
        for k = 1:runSize
            startInd=((k-1)*numConds)+1;
            endInd=k*numConds;
            actMat(:,startInd:endInd)=vals(tempSLinds{k,1},:);
        end
        if useGPU==1
            RSMs(:,selectInd)= CreateBatchedRSMsGPU(actMat,numConds,runSize,RSMsize); 
        else
            RSMs(:,selectInd)= CreateBatchedRSMs(actMat,numConds,runSize,RSMsize);
        end
    else
        batch_StartInd=[1:MaxBatchSize:runSize]';
        batch_EndInd=[batch_StartInd(2:end,1)-1;runSize];            
        numBatches=length(batch_StartInd);
        for j = 1:numBatches
            batch_selectInd=selectInd(batch_StartInd(j):batch_EndInd(j));
            batch_runSize=length(batch_selectInd);
            tempSLinds=SLinds(batch_selectInd,:);
            actMat=zeros(i,batch_runSize*10,'single');
            for k = 1:batch_runSize
                startInd=((k-1)*numConds)+1;
                endInd=k*numConds;
                actMat(:,startInd:endInd)=vals(tempSLinds{k,1},:);
            end   
            if useGPU==1
                RSMs(:,batch_selectInd)= CreateBatchedRSMsGPU(actMat,numConds,batch_runSize,RSMsize); 
            else
                RSMs(:,batch_selectInd)= CreateBatchedRSMs(actMat,numConds,batch_runSize,RSMsize); 
            end
        end
    end
end
RSMs=gather(RSMs);
end

function RSMs = CreateBatchedRSMsGPU(actMat,numConds,runSize,RSMsize)
    SelectMat=repmat({ones(numConds)},[1,runSize]);
    SelectMat=blkdiag(SelectMat{:});
    SelectMat=triu(SelectMat,1)==1;
    actMat=corrcoef(gpuArray(actMat));
    RSMs=reshape(actMat(SelectMat),[RSMsize,runSize]);
end
function RSMs = CreateBatchedRSMs(actMat,numConds,runSize,RSMsize)
    SelectMat=repmat({ones(numConds)},[1,runSize]);
    SelectMat=blkdiag(SelectMat{:});
    SelectMat=triu(SelectMat,1)==1;
    actMat=corrcoef(actMat);
    RSMs=reshape(actMat(SelectMat),[RSMsize,runSize]);
end