function [OutVars,GroupData] = FullRCA(Data,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vertIn = VariableSetter( 'vertIn',1,varargin);
MakeSymmetric = VariableSetter( 'MakeSymmetric',1,varargin);
RSMGroup = VariableSetter( 'RSMGroup','All',varargin);
DownSample = VariableSetter( 'DownSample',0,varargin);
DownSampleN = VariableSetter( 'DownSampleN',[],varargin);
DownSampleReps = VariableSetter( 'DownSampleReps',100,varargin);
RunPermute = VariableSetter( 'RunPermute',0,varargin);
PermuteReps = VariableSetter( 'PermuteReps',1000,varargin);
OutputPermDists = VariableSetter( 'OutputPermDists',0,varargin);
zNorm = VariableSetter( 'zNorm',0,varargin);
if vertIn==1
    numSS=size(Data,3);
else
    numSS=size(Data,4);
end
if DownSample==1
    if numSS<=DownSampleN
        disp('Down sample size is not less than overall N: skipping')
        DownSample=0;
    end
end

if strcmpi(RSMGroup,'W')
    [Data,~] = RSM2RSMGroups(Data,vertIn);
elseif strcmpi(RSMGroup,'X')
    [~,Data] = RSM2RSMGroups(Data,vertIn);
end

[~,OutVars.RF,OutVars.RCMat,GroupData]=FastRCA( Data,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );

if DownSample == 1
    DownSample_RF=nan(size(OutVars.RF,1),size(OutVars.RF,2),DownSampleReps,'single');
  	DownSample_RC=nan(size(OutVars.RCMat,1),size(OutVars.RCMat,2),size(OutVars.RCMat,3),DownSampleReps,'single'); 
    
    parfor subRep=1:DownSampleReps
        tempDownSample_RF=nan(size(OutVars.RF,1),size(OutVars.RF,2));
        tempDownSample_RC=nan(size(OutVars.RCMat,1),size(OutVars.RCMat,2),size(OutVars.RCMat,3));
        subSampleInd=randperm(numSS,DownSampleN);
        if vertIn==1
            subData=Data(:,:,subSampleInd);
        else
            subData=Data(:,:,:,subSampleInd);
        end
        [~,tempDownSample_RF(:,subSampleInd),tempDownSample_RC(:,:,subSampleInd)]=FastRCA( subData,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
        DownSample_RF(:,:,subRep)=tempDownSample_RF;
        DownSample_RC(:,:,:,subRep)=tempDownSample_RC;
    end
    OutVars.RF=tanh(nanmean(atanh(DownSample_RF),3)); 
    OutVars.RCMat=tanh(nanmean(atanh(DownSample_RC),4));  
end
if RunPermute == 1
    [ ~,RealRF,~,RealRC ] = FastRCA( Data,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
    RealRF=mean(atanh(RealRF),2)';
    RealRC=mean(atanh(RealRC),2)';
    [RF_Dist,RC_Dist]=PermuteRCA(Data,PermuteReps,vertIn,MakeSymmetric,zNorm);
    if OutputPermDists==1
        OutVars.RF_Dist=RF_Dist;
        OutVars.RC_Dist=RC_Dist;
    else
        OutVars.RF_Dist=[];
        OutVars.RC_Dist=[];
    end
    RF_Dist=mean(atanh(RF_Dist),3);
    RC_Dist=abs(mean(atanh(RC_Dist),3));
    [RF_perctiles] = reverseprctile(RF_Dist,RealRF);
    RF_perctiles=1-RF_perctiles;    
    [RF_permZs]=p2z(RF_perctiles,2);    
    [RC_perctiles] = reverseprctile(RC_Dist,abs(RealRC));
    RC_perctiles=1-RC_perctiles;
    RC_permZs=p2z(RC_perctiles,2).*sign(RealRC);
    OutVars.RF_permP=RF_perctiles;
    OutVars.RF_permZ=RF_permZs; 
    OutVars.RCMat_permP = vertRSM2SymRSM( RC_perctiles(:));
    OutVars.RCMat_permZ = vertRSM2SymRSM( RC_permZs(:));  
    
    RF_Dist=max(RF_Dist,[],2);
    RC_Dist=max(RC_Dist,[],2);
    [RF_perctiles] = reverseprctile(RF_Dist,RealRF);
    RF_perctiles=1-RF_perctiles;    
    [RF_permZs]=p2z(RF_perctiles,2);    
    [RC_perctiles] = reverseprctile(RC_Dist,abs(RealRC));
    RC_perctiles=1-RC_perctiles;
    RC_permZs=p2z(RC_perctiles,2).*sign(RealRC);
    OutVars.RF_permP_Corrected=RF_perctiles;
    OutVars.RF_permZ_Corrected=RF_permZs; 
    OutVars.RCMat_permP_Corrected = vertRSM2SymRSM( RC_perctiles(:));
    OutVars.RCMat_permZ_Corrected = vertRSM2SymRSM( RC_permZs(:));  
    
else
    OutVars.RF_permP=[];
    OutVars.RF_permZ=[]; 
    OutVars.RCMat_permP = [];
    OutVars.RCMat_permZ = [];    
    OutVars.RF_permP_Corrected=[];
    OutVars.RF_permZ_Corrected=[]; 
    OutVars.RCMat_permP_Corrected = [];
    OutVars.RCMat_permZ_Corrected = [];     
end
end

function [RF_Dist,RC_Dist]=PermuteRCA(Data,numReps,vertIn,MakeSymmetric,zNorm)

if vertIn == 1  
    numSS = size(Data,3);
    numROIs = size(Data,2);
    tData=[];
    for ssNum=1:numSS
        [ SymRSMs ] = vertRSM2SymRSM( Data(:,:,ssNum ));
        tData=cat(4,tData,SymRSMs);
    end
    Data=tData;
else
    numSS = size(Data,4);
    numROIs = size(Data,3);
end
numCxns=(numROIs^2-numROIs)/2;
RF_Dist=zeros(numROIs,numSS,numReps,'single');
RC_Dist=zeros(numCxns,numSS,numReps,'single');
parfor rep = 1:numReps
    permData=Data*0;
    for ssNum=1:numSS
        [permData(:,:,:,ssNum)]=modelscramble(Data(:,:,:,ssNum ),'VH_once');
    end
    [ ~,RF_Dist(:,:,rep),~,RC_Dist(:,:,rep) ] = FastRCA( permData,'vertIn',0,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
end
RF_Dist=permute(RF_Dist,[3,1,2]);
RC_Dist=permute(RC_Dist,[3,1,2]);
end