function [ OutVars,RF,RCMat,RCVert,GroupData ] = FastRCA( Data,varargin )
%Written by David Rothlein
% Compute cross-validated representational fidelity and representational connectivity values 
vertIn = VariableSetter( 'vertIn',1,varargin);
MakeSymmetric = VariableSetter( 'MakeSymmetric',1,varargin);
corrType = VariableSetter( 'corrType','Pearson',varargin);
CrossValType = VariableSetter( 'CrossValType','LeaveOneOut',varargin);
zNorm = VariableSetter( 'zNorm',0,varargin);
zNorm_Post = VariableSetter( 'zNorm_Post',0,varargin);
if strcmpi(CrossValType,'LeaveOneOut')
    if vertIn==1
        RSMvert=Data;
        numROIs=size(RSMvert,2);
        numSS=size(RSMvert,3); 
    else    
        [RSMvert,RSMs] = RSMformatConvert( Data );
        numROIs=size(RSMs,3);
        numSS=size(RSMs,4);        
    end
    if zNorm==1
        for i = 1:numSS
            RSMvert(:,:,i)=zscore(RSMvert(:,:,i));
        end
    end
    GroupData=repmat(nansum(RSMvert,3)/(numSS-1),[1,1,numSS])-(RSMvert/(numSS-1));
    OutVars.RF=zeros(numROIs,numSS,'single');
    OutVars.RCMat=zeros(numROIs,numROIs,numSS,'single');
    VertLength=(numROIs^2-numROIs)/2;
    if MakeSymmetric==1
        OutVars.RCVert=zeros(VertLength,numSS,'single');
    else
        OutVars.RCVert=zeros(VertLength,numSS,2,'single');
    end

    RFInd= eye(numROIs);
    TriuInd=triu(ones(numROIs),1);
    if strcmpi(corrType,'Pearson') || strcmpi(corrType,'Spearman')
        for ss=1:numSS
            FullMat=corr(RSMvert(:,:,ss),GroupData(:,:,ss),'type',corrType);
            OutVars.RF(:,ss)=FullMat(RFInd==1);    
            FullMat(RFInd==1)=nan;
            if MakeSymmetric==1
                FullMat=(FullMat+permute(FullMat,[2,1]))/2;
                if zNorm_Post == 1
                    FullMat=vertRSM2SymRSM(nan_zscore(mat2uppertriuvectormat(FullMat)),[],1);
                end
                OutVars.RCMat(:,:,ss)=FullMat;
                OutVars.RCVert(:,ss)=FullMat(TriuInd==1);
            else
                if zNorm_Post == 1
                    FullMat=nan_zscore(FullMat,'pooled');
                end
                OutVars.RCMat(:,:,ss)=FullMat;
                OutVars.RCVert(:,ss,1)=FullMat(TriuInd==1);
                FullMat=permute(FullMat,[2,1]);
                OutVars.RCVert(:,ss,2)=FullMat(TriuInd==1);
            end      
        end
    end
    if zNorm_Post == 1
        OutVars.RF=nan_zscore(OutVars.RF);
    end
    RF=OutVars.RF;
    RCMat=OutVars.RCMat;
    RCVert=OutVars.RCVert;
end

end