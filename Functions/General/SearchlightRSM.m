function [ RSMs,meanRSMs,rsm_mask,SLinds,numVoxels,MaskCoords] = SearchlightRSM( vals,Mask,SLrad,SLThresh,Shape,MaxBatchSize,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%BatchSize of 5 = fastest
%vals= activation patterns (Voxel,Condition)
[useGPU] = VariableSetter('useGPU',0,varargin);
[SLinds] = VariableSetter('SLinds',[],varargin);
[numVoxels] = VariableSetter('numVoxels',[],varargin);
[MaskCoords] = VariableSetter('MaskCoords',[],varargin);
ComputeMeanRSM = VariableSetter('ComputeMeanRSM',0,varargin);

numConds=size(vals,2);
RSMsize=(numConds^2-numConds)/2;
%Generate searchlight indicies from mask
if isempty(SLinds)
    [SLinds,numVoxels,MaskCoords] = Mask2Searchlight(Mask,SLrad,SLThresh,Shape);
    CleanInd = numVoxels < SLThresh;
    SLinds(CleanInd,:)=[];
    MaskCoords(CleanInd,:)=[];
    numVoxels(CleanInd,:)=[];
end
CleanInd = numVoxels < SLThresh;
SLinds(CleanInd,:)=[];
MaskCoords(CleanInd,:)=[];
numVoxels(CleanInd,:)=[];
numSLs=length(SLinds);
SLsizes=unique(numVoxels);

if ComputeMeanRSM == 1
    meanRSMs=zeros(RSMsize,numSLs,'single');
else
    meanRSMs=[];
end
if useGPU==1
    gpuDevice;
    RSMs=zeros(RSMsize,numSLs,'gpuArray');
else
    RSMs=zeros(RSMsize,numSLs,'single');
end
rsm_mask=coords2mat(MaskCoords,Mask,ones(numSLs,1));
vertMask=Mask(:);
vertInd=vertMask*0;
skipSLCount=0;
for i = SLsizes'
    selectInd=find(numVoxels==i);
    runSize=length(selectInd);       
    if runSize<=MaxBatchSize
        tempSLinds=SLinds(selectInd,:);
        actMat=zeros(i,runSize*numConds,'single');
        for k = 1:runSize
            startInd=((k-1)*numConds)+1;
            endInd=k*numConds;
            TempVertInd=vertInd;
            TempVertInd(tempSLinds{k,1},1)=1;
            TempVertInd=TempVertInd(vertMask~=0,:)~=0;            
            actMat(:,startInd:endInd)=vals(TempVertInd,:);
        end
        if useGPU==1
            RSMs(:,selectInd)= CreateBatchedRSMsGPU(actMat,numConds,runSize,RSMsize); 
        else
            RSMs(:,selectInd)= CreateBatchedRSMs(actMat,numConds,runSize,RSMsize);
        end
        if ComputeMeanRSM == 1
            meanRSMs(:,selectInd)= CreateBatchedMeanRSMs(nanmean(actMat,1),numConds,runSize,RSMsize);
        end
    else
        batch_StartInd=[1:MaxBatchSize:runSize]';
        batch_EndInd=[batch_StartInd(2:end,1)-1;runSize];            
        numBatches=length(batch_StartInd);
        for j = 1:numBatches
            batch_selectInd=selectInd(batch_StartInd(j):batch_EndInd(j));
            batch_runSize=length(batch_selectInd);
            tempSLinds=SLinds(batch_selectInd,:);
            actMat=zeros(i,batch_runSize*numConds,'single');
            removeInd=zeros(batch_runSize,1);
            removeInd_ActMat=zeros(1,batch_runSize*numConds);
            for k = 1:batch_runSize
                startInd=((k-1)*numConds)+1;
                endInd=k*numConds;
                TempVertInd=vertInd;
                TempVertInd(tempSLinds{k,1},1)=1;
                TempVertInd=TempVertInd(vertMask~=0,:)~=0;
                if sum(TempVertInd)~=i
                    removeInd(k,1)=1;
                    removeInd_ActMat(1,startInd:endInd)=ones(1,numConds);
                    skipSLCount=skipSLCount+1;
                    continue
                else
                    actMat(:,startInd:endInd)=vals(TempVertInd,:);
                end
            end   
            batch_selectInd(removeInd==1)=[];
            batch_runSize=length(batch_selectInd);
            actMat(:,removeInd_ActMat==1)=[];
            if useGPU==1
                RSMs(:,batch_selectInd)= CreateBatchedRSMsGPU(actMat,numConds,batch_runSize,RSMsize); 
            else
                RSMs(:,batch_selectInd)= CreateBatchedRSMs(actMat,numConds,batch_runSize,RSMsize); 
            end
            if ComputeMeanRSM == 1
                meanRSMs(:,batch_selectInd)= CreateBatchedMeanRSMs(nanmean(actMat,1),numConds,batch_runSize,RSMsize);
            end            
        end
    end
end
if useGPU==1
    RSMs=gather(RSMs);
end
%disp(['SLs Skipped = ',num2str(skipSLCount)]);
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
function meanRSMs = CreateBatchedMeanRSMs(actMat,numConds,runSize,RSMsize)
    SelectMat=repmat({ones(numConds)},[1,runSize]);
    SelectMat=blkdiag(SelectMat{:});
    SelectMat=triu(SelectMat,1)==1;
    [actMat] = vector2diffdistmat( actMat );
    meanRSMs=reshape(actMat(SelectMat),[RSMsize,runSize]);
end
function [diffpropmat] = vector2diffdistmat( v )
    diffmat=abs(bsxfun(@minus,v,v'));
    maxv=max(diffmat(:));
    diffpropmat=1-(diffmat/maxv);
end
