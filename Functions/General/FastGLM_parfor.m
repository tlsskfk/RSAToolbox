function [Condition_tMat,DesignMatrix,Condition_bMat] = FastGLM_parfor(Y,X,contrasts,batchSize,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numY=size(Y,2);
numX=size(X,2);
[EmptyRegressors] = VariableSetter('EmptyRegressors',zeros(1,numX),varargin);
[ResampleToFit] = VariableSetter('ResampleToFit','Y',varargin);
[UseStat] = VariableSetter('UseStat','T',varargin);
UseRegressors=EmptyRegressors==0;
ResampDur=size(X,1);
skipContrast=0;
if ~isempty(contrasts) && sum(EmptyRegressors(:))>0
    UseContrast=sum(abs(contrasts(:,EmptyRegressors(1,[1:size(contrasts,2)])==1)),2)==0;
    tMat_UseInd=[UseRegressors,single(UseContrast)'];
    if sum(single(UseContrast(:)))==0
        skipContrast=1;
    end
elseif ~isempty(contrasts)
    UseContrast=zeros(size(contrasts,1),1)==0;
    tMat_UseInd=[UseRegressors,single(UseContrast)'];
else
    tMat_UseInd=UseRegressors;
end
tMat_UseInd=tMat_UseInd==1;
if size(Y,1)~=size(X,1)
    if strcmpi(ResampleToFit,'Y')
        X=imresize(X,[size(Y,1),size(X,2)]);
    else
        Y=imresize(Y,[size(X,1),size(Y,2)]);
    end
end
DesignMatrix=X(:,UseRegressors==1);
if batchSize<numY
    batch_StartInd=[1:batchSize:numY];
    batch_EndInd=[batch_StartInd(1,2:end)-1,numY];
    numBatches=length(batch_StartInd);
    tempCondition_tMat=cell(numBatches,1);
    tempCondition_bMat=cell(numBatches,1);
    batch_Ys=cell(numBatches,1);
    for batch=1:numBatches
        batch_Ys{batch,1}=Y(:,[batch_StartInd(1,batch):batch_EndInd(1,batch)]);
    end
    if ~isempty(contrasts)
        numContrasts=size(contrasts,1);
        if size(contrasts,2)<numX
            contrasts=[contrasts,zeros(numContrasts,numX-size(contrasts,2))]; 
        end    
        Condition_tMat=single(zeros(numY,numX+numContrasts));
        Condition_bMat=single(zeros(numY,numX+numContrasts));
        parfor batch=1:numBatches
            batch_Y=batch_Ys{batch,1};
            if skipContrast==0
                [T,B] = FastOLSRegressContrast(batch_Y,X(:,UseRegressors==1),contrasts(UseContrast==1,UseRegressors==1));
            else
                [T,B] = FastOLSRegress(batch_Y,X(:,UseRegressors==1));
            end
            tempCondition_bMat{batch,1}=B;
            tempCondition_tMat{batch,1}=T;  
            warning('off','all'); 
        end
    else  
        Condition_bMat=single(zeros(numY,numX));
        Condition_tMat=single(zeros(numY,numX));
        parfor batch=1:numBatches
            batch_Y=batch_Ys{batch,1};
            [T,B] = FastOLSRegress(batch_Y,X(:,UseRegressors==1));
            tempCondition_bMat{batch,1}=B;
            tempCondition_tMat{batch,1}=T;  
            warning('off','all');
        end
    end
    for batch=1:numBatches
        Condition_tMat([batch_StartInd(1,batch):batch_EndInd(1,batch)],tMat_UseInd)=tempCondition_tMat{batch,1};
        Condition_bMat([batch_StartInd(1,batch):batch_EndInd(1,batch)],tMat_UseInd)=tempCondition_bMat{batch,1};
    end    
else
    if ~isempty(contrasts)
        numContrasts=size(contrasts,1);
        Condition_tMat=single(zeros(numY,numX+numContrasts));
        Condition_bMat=single(zeros(numY,numX+numContrasts));
        if size(contrasts,2)<numX
            contrasts=[contrasts,zeros(numContrasts,numX-size(contrasts,2))]; 
        end   
        if skipContrast==0
            [T,B] = FastOLSRegressContrast(Y,X(:,UseRegressors==1),contrasts(UseContrast==1,UseRegressors==1));
        else
            [T,B] = FastOLSRegress(Y,X(:,UseRegressors==1));
        end
        warning('off','all');
        Condition_bMat(:,tMat_UseInd)=B;
        Condition_tMat(:,tMat_UseInd)=T;
        
    else  
        Condition_tMat=single(zeros(numY,numX));
        Condition_bMat=single(zeros(numY,numX));
        [T,B] = FastOLSRegress(Y,X(:,UseRegressors==1));          
        Condition_bMat(:,tMat_UseInd)=B;
        Condition_tMat(:,tMat_UseInd)=T;
                
        warning('off','all');
    end
end

