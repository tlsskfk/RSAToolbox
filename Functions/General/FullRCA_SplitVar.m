function [OutVars1,OutVars2,OutVarsDiff] = FullRCA_SplitVar(Data1,Data2,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vertIn = VariableSetter( 'vertIn',1,varargin);
MakeSymmetric = VariableSetter( 'MakeSymmetric',1,varargin);
DownSample = VariableSetter( 'DownSample',0,varargin);
DownSampleN = VariableSetter( 'DownSampleN',[],varargin);
DownSampleReps = VariableSetter( 'DownSampleReps',100,varargin);
RunPermute = VariableSetter( 'RunPermute',0,varargin);
PermuteReps = VariableSetter( 'PermuteReps',1000,varargin);
OutputPermDists = VariableSetter( 'OutputPermDists',0,varargin);
zNorm = VariableSetter( 'zNorm',0,varargin);
if vertIn==1
    numSS=size(Data1,3);
else
    numSS=size(Data1,4);
end
if DownSample==1
    if numSS<=DownSampleN
        disp('Down sample size is not less than overall N: skipping')
        DownSample=0;
    end
end

[~,OutVars1.RF,OutVars1.RCMat]=FastRCA( Data1,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm);
[~,OutVars2.RF,OutVars2.RCMat]=FastRCA( Data2,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm);
OutVarsDiff.RF=tanh(atanh(OutVars2.RF)-atanh(OutVars1.RF));
OutVarsDiff.RCMat=tanh(atanh(OutVars2.RCMat)-atanh(OutVars1.RCMat));

if DownSample == 1
    DownSample_RF1=nan(size(OutVars1.RF,1),size(OutVars1.RF,2),DownSampleReps,'single');
  	DownSample_RC1=nan(size(OutVars1.RCMat,1),size(OutVars1.RCMat,2),size(OutVars1.RCMat,3),DownSampleReps,'single'); 
    DownSample_RF2=nan(size(OutVars2.RF,1),size(OutVars2.RF,2),DownSampleReps,'single');
  	DownSample_RC2=nan(size(OutVars2.RCMat,1),size(OutVars2.RCMat,2),size(OutVars2.RCMat,3),DownSampleReps,'single');     
    parfor subRep=1:DownSampleReps
        tempDownSample_RF1=nan(size(OutVars1.RF,1),size(OutVars1.RF,2));
        tempDownSample_RC1=nan(size(OutVars1.RCMat,1),size(OutVars1.RCMat,2),size(OutVars1.RCMat,3));
        tempDownSample_RF2=nan(size(OutVars2.RF,1),size(OutVars2.RF,2));
        tempDownSample_RC2=nan(size(OutVars2.RCMat,1),size(OutVars2.RCMat,2),size(OutVars2.RCMat,3));        
        subSampleInd=randperm(numSS,DownSampleN);
        if vertIn==1
            subData1=Data1(:,:,subSampleInd);
            subData2=Data2(:,:,subSampleInd);
        else
            subData1=Data1(:,:,:,subSampleInd);
            subData2=Data2(:,:,:,subSampleInd);
        end
        [~,tempDownSample_RF1(:,subSampleInd),tempDownSample_RC1(:,:,subSampleInd)]=FastRCA( subData1,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
        DownSample_RF1(:,:,subRep)=tempDownSample_RF1;
        DownSample_RC1(:,:,:,subRep)=tempDownSample_RC1;
        [~,tempDownSample_RF2(:,subSampleInd),tempDownSample_RC2(:,:,subSampleInd)]=FastRCA( subData2,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
        DownSample_RF2(:,:,subRep)=tempDownSample_RF2;
        DownSample_RC2(:,:,:,subRep)=tempDownSample_RC2;        
    end
    OutVars1.RF=tanh(nanmean(atanh(DownSample_RF1),3)); 
    OutVars1.RCMat=tanh(nanmean(atanh(DownSample_RC1),4));  
    OutVars2.RF=tanh(nanmean(atanh(DownSample_RF2),3)); 
    OutVars2.RCMat=tanh(nanmean(atanh(DownSample_RC2),4)); 
    OutVarsDiff.RF=tanh(nanmean(atanh(DownSample_RF2)-atanh(DownSample_RF1),3)); 
    OutVarsDiff.RCMat=tanh(nanmean(atanh(DownSample_RC2)-atanh(DownSample_RC1),4));     
end

if RunPermute == 1
    [ ~,RealRF1,~,RealRC1 ] = FastRCA( Data1,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
    [ ~,RealRF2,~,RealRC2 ] = FastRCA( Data2,'vertIn',vertIn,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
    RealRF1=tanh(mean(atanh(RealRF1),2))';
    RealRC1=tanh(mean(atanh(RealRC1),2))';
    RealRF2=tanh(mean(atanh(RealRF2),2))';
    RealRC2=tanh(mean(atanh(RealRC2),2))'; 
    RealRFDiff=tanh(atanh(RealRF2)-atanh(RealRF1));
    RealRCDiff=tanh(atanh(RealRC2)-atanh(RealRC1));  
    
    [RF_Dist1,RC_Dist1,rperms]=PermuteRCA(Data1,PermuteReps,vertIn,MakeSymmetric,[],zNorm);
    [RF_Dist2,RC_Dist2]=PermuteRCA(Data2,PermuteReps,vertIn,MakeSymmetric,rperms,zNorm);    
    RF_DistDiff=tanh(atanh(RF_Dist2)-atanh(RF_Dist1));
    RC_DistDiff=tanh(atanh(RC_Dist2)-atanh(RC_Dist1));
    if OutputPermDists==1
        OutVars1.RF_Dist=RF_Dist1;
        OutVars1.RC_Dist=RC_Dist1;
    else
        OutVars1.RF_Dist=[];
        OutVars1.RC_Dist=[];
    end
    if OutputPermDists==1
        OutVars2.RF_Dist=RF_Dist2;
        OutVars2.RC_Dist=RC_Dist2;
    else
        OutVars2.RF_Dist=[];
        OutVars2.RC_Dist=[];
    end
    if OutputPermDists==1
        OutVarsDiff.RF_Dist=RF_DistDiff;
        OutVarsDiff.RC_Dist=RC_DistDiff;
    else
        OutVarsDiff.RF_Dist=[];
        OutVarsDiff.RC_Dist=[];
    end   
    
    RF_Dist1=tanh(mean(atanh(RF_Dist1),3));
    RC_Dist1=tanh(abs(mean(atanh(RC_Dist1),3)));
    [RF_perctiles1] = reverseprctile(RF_Dist1,RealRF1);
    RF_perctiles1=1-RF_perctiles1;    
    [RF_permZs1]=p2z(RF_perctiles1,2);    
    [RC_perctiles1] = reverseprctile(RC_Dist1,abs(RealRC1));
    RC_perctiles1=1-RC_perctiles1;
    RC_permZs1=p2z(RC_perctiles1,2).*sign(RealRC1);
    OutVars1.RF_permP=RF_perctiles1;
    OutVars1.RF_permZ=RF_permZs1; 
    OutVars1.RCMat_permP = vertRSM2SymRSM( RC_perctiles1(:));
    OutVars1.RCMat_permZ = vertRSM2SymRSM( RC_permZs1(:));  

    RF_Dist2=tanh(mean(atanh(RF_Dist2),3));
    RC_Dist2=tanh(abs(mean(atanh(RC_Dist2),3)));
    [RF_perctiles2] = reverseprctile(RF_Dist2,RealRF2);
    RF_perctiles2=1-RF_perctiles2;    
    [RF_permZs2]=p2z(RF_perctiles2,2);    
    [RC_perctiles2] = reverseprctile(RC_Dist2,abs(RealRC2));
    RC_perctiles2=1-RC_perctiles2;
    RC_permZs2=p2z(RC_perctiles2,2).*sign(RealRC2);
    OutVars2.RF_permP=RF_perctiles2;
    OutVars2.RF_permZ=RF_permZs2; 
    OutVars2.RCMat_permP = vertRSM2SymRSM( RC_perctiles2(:));
    OutVars2.RCMat_permZ = vertRSM2SymRSM( RC_permZs2(:)); 
    
    RF_DistDiff=tanh(abs(mean(atanh(RF_DistDiff),3)));
    RC_DistDiff=tanh(abs(mean(atanh(RC_DistDiff),3)));
    [RF_perctilesDiff] = reverseprctile(RF_DistDiff,abs(RealRFDiff));
    RF_perctilesDiff=1-RF_perctilesDiff;    
    [RF_permZsDiff]=p2z(RF_perctilesDiff,2).*sign(RealRFDiff);    
    [RC_perctilesDiff] = reverseprctile(RC_DistDiff,abs(RealRCDiff));
    RC_perctilesDiff=1-RC_perctilesDiff;
    RC_permZsDiff=p2z(RC_perctilesDiff,2).*sign(RealRCDiff);
    OutVarsDiff.RF_permP=RF_perctilesDiff;
    OutVarsDiff.RF_permZ=RF_permZsDiff; 
    OutVarsDiff.RCMat_permP = vertRSM2SymRSM( RC_perctilesDiff(:));
    OutVarsDiff.RCMat_permZ = vertRSM2SymRSM( RC_permZsDiff(:));     
    
    RF_Dist1=max(RF_Dist1,[],2);
    RC_Dist1=max(RC_Dist1,[],2);
    [RF_perctiles1] = reverseprctile(RF_Dist1,RealRF1);
    RF_perctiles1=1-RF_perctiles1;    
    [RF_permZs1]=p2z(RF_perctiles1,2);    
    [RC_perctiles1] = reverseprctile(RC_Dist1,abs(RealRC1));
    RC_perctiles1=1-RC_perctiles1;
    RC_permZs1=p2z(RC_perctiles1,2).*sign(RealRC1);
    OutVars1.RF_permP_Corrected=RF_perctiles1;
    OutVars1.RF_permZ_Corrected=RF_permZs1; 
    OutVars1.RCMat_permP_Corrected = vertRSM2SymRSM( RC_perctiles1(:));
    OutVars1.RCMat_permZ_Corrected = vertRSM2SymRSM( RC_permZs1(:));  
    
    RF_Dist2=max(RF_Dist2,[],2);
    RC_Dist2=max(RC_Dist2,[],2);
    [RF_perctiles2] = reverseprctile(RF_Dist2,RealRF2);
    RF_perctiles2=1-RF_perctiles2;    
    [RF_permZs2]=p2z(RF_perctiles2,2);    
    [RC_perctiles2] = reverseprctile(RC_Dist2,abs(RealRC2));
    RC_perctiles2=1-RC_perctiles2;
    RC_permZs2=p2z(RC_perctiles2,2).*sign(RealRC2);
    OutVars2.RF_permP_Corrected=RF_perctiles2;
    OutVars2.RF_permZ_Corrected=RF_permZs2; 
    OutVars2.RCMat_permP_Corrected = vertRSM2SymRSM( RC_perctiles2(:));
    OutVars2.RCMat_permZ_Corrected = vertRSM2SymRSM( RC_permZs2(:)); 
    
    RF_DistDiff=max(RF_DistDiff,[],2);
    RC_DistDiff=max(RC_DistDiff,[],2);
    [RF_perctilesDiff] = reverseprctile(RF_DistDiff,abs(RealRFDiff));
    RF_perctilesDiff=1-RF_perctilesDiff;    
    [RF_permZsDiff]=p2z(RF_perctilesDiff,2).*sign(RealRFDiff);    
    [RC_perctilesDiff] = reverseprctile(RC_DistDiff,abs(RealRCDiff));
    RC_perctilesDiff=1-RC_perctilesDiff;
    RC_permZsDiff=p2z(RC_perctilesDiff,2).*sign(RealRCDiff);
    OutVarsDiff.RF_permP_Corrected=RF_perctilesDiff;
    OutVarsDiff.RF_permZ_Corrected=RF_permZsDiff; 
    OutVarsDiff.RCMat_permP_Corrected = vertRSM2SymRSM( RC_perctilesDiff(:));
    OutVarsDiff.RCMat_permZ_Corrected = vertRSM2SymRSM( RC_permZsDiff(:));     
    
else
    OutVars1.RF_permP=[];
    OutVars1.RF_permZ=[]; 
    OutVars1.RCMat_permP = [];
    OutVars1.RCMat_permZ = [];    
    OutVars1.RF_permP_Corrected=[];
    OutVars1.RF_permZ_Corrected=[]; 
    OutVars1.RCMat_permP_Corrected = [];
    OutVars1.RCMat_permZ_Corrected = [];
    
    OutVars2.RF_permP=[];
    OutVars2.RF_permZ=[]; 
    OutVars2.RCMat_permP = [];
    OutVars2.RCMat_permZ = [];    
    OutVars2.RF_permP_Corrected=[];
    OutVars2.RF_permZ_Corrected=[]; 
    OutVars2.RCMat_permP_Corrected = [];
    OutVars2.RCMat_permZ_Corrected = [];  
    
    OutVarsDiff.RF_permP=[];
    OutVarsDiff.RF_permZ=[]; 
    OutVarsDiff.RCMat_permP = [];
    OutVarsDiff.RCMat_permZ = [];    
    OutVarsDiff.RF_permP_Corrected=[];
    OutVarsDiff.RF_permZ_Corrected=[]; 
    OutVarsDiff.RCMat_permP_Corrected = [];
    OutVarsDiff.RCMat_permZ_Corrected = [];      
end
end

function [RF_Dist,RC_Dist,rperms]=PermuteRCA(Data,numReps,vertIn,MakeSymmetric,rperms,zNorm)

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
        [permData(:,:,:,ssNum)]=modelscramble(Data(:,:,:,ssNum ),'VH_once',rperms);
    end
    [ ~,RF_Dist(:,:,rep),~,RC_Dist(:,:,rep) ] = FastRCA( permData,'vertIn',0,'MakeSymmetric',MakeSymmetric,'zNorm',zNorm );
end
RF_Dist=permute(RF_Dist,[3,1,2]);
RC_Dist=permute(RC_Dist,[3,1,2]);
end