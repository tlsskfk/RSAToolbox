function [RC,RF,IF,RCs,RCx,RFs,RFx,rsm_mask] = GetParcellationSplitRCA(Events,ConfoundEvents,TimeCourses,ConfoundTCs,boldTCs,boldMasks,ConditionNames,Parcels,TrialNum,NumReps,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%% Optional input for GetParcellationSplitRCA
[SplitAcrossRun] = VariableSetter('SplitAcrossRun',0,varargin);
%%% Optional input for RSMParcellationRSM
[csfMask] = VariableSetter('csfMask',[],varargin);
[wmMask] = VariableSetter('wmMask',[],varargin);
[Compute_MeanAct_RSM] = VariableSetter('Compute_MeanAct_RSM',0,varargin);
[gmMask] = VariableSetter('gmMask',[],varargin);
[SimType] = VariableSetter('SimType','corrcoef',varargin);
[OutputConfoundRSMs] = VariableSetter('OutputConfoundRSMs',0,varargin);

%%% Optional input for GetActivationPatterns
[parGLM] = VariableSetter('parGLM',1,varargin);
[Compute_CondReg_CorrMat] = VariableSetter('Compute_CondReg_CorrMat',0,varargin);
[ContrastNames] = VariableSetter('ContrastNames',[],varargin);
%%% Optional input for GLM_RegressorPrep
[Normalize] = VariableSetter('Normalize','zscore',varargin);
[ConfoundTCsBySubj] = VariableSetter('ConfoundTCsBySubj',{[]},varargin);
[ResampleSizes] = VariableSetter('ResampleSizes',[],varargin);
[AM2]=VariableSetter('AM2',1,varargin);
[AM2_Events]=VariableSetter('AM2_Events',{[]},varargin);
TimecourseShift = VariableSetter('TimecourseShift',0,varargin);

%%% Optional input for BoldTC
[ParcellationMask] = VariableSetter('ParcellationMask',[],varargin);
[ResampleSizesBoldTC] = VariableSetter('ResampleSizesBoldTC',[],varargin);
[NormWithinRun] = VariableSetter('NormWithinRun','',varargin);
[NormAcrossRun] = VariableSetter('NormAcrossRun','',varargin);

%%% Optional input for FastGLM
[Contrasts] = VariableSetter('Contrasts',[],varargin);
[BatchSize] = VariableSetter('BatchSize',1000,varargin);
[ResampleToFit] = VariableSetter('ResampleToFit','Y',varargin);
[UseStat] = VariableSetter('UseStat','T',varargin); %T of B

numRuns=size(Events,1);
numParcels=size(Parcels,1);

Temp_RC=cell(NumReps,1);
Temp_RF=cell(NumReps,1);
Temp_IF=cell(NumReps,1);
Temp_RCs=cell(NumReps,1);
Temp_RFs=cell(NumReps,1);
Temp_RCx=cell(NumReps,1);
Temp_RFx=cell(NumReps,1);

All_RC=cell(numParcels,1);
All_RF=cell(numParcels,1);
All_IF=cell(numParcels,1);
All_RCs=cell(numParcels,1);
All_RFs=cell(numParcels,1);
All_RCx=cell(numParcels,1);
All_RFx=cell(numParcels,1);
if SplitAcrossRun==1
    [EventsMat,CatInd] = Cell2CatMat(Events(:,1));
    TrialMat=[];
    for run=1:numRuns
        if run > 1
            MaxVal=max(TrialMat(:));
            TrialMat=[TrialMat;TrialNum{run,1}+MaxVal];
        else
            TrialMat=TrialNum{run,1};
        end
    end
else
    EventsMat=[];
    TrialMat=[];
    CatInd=[];
end
rsm_mask=cell(numParcels,NumReps);
parfor rep = 1:NumReps
    Use_Events=cell(numRuns,2);
    if SplitAcrossRun==1
        [SplitConditions,SplitConditionNames,~] = ComputeSplitConditions(EventsMat,ConditionNames,'SplitVar',rand(size(EventsMat,1),1),'TrialNum',TrialMat,'RandomSplit',1);
        [SplitConditions] = CatMat2Cell(SplitConditions,CatInd);        
        for run=1:numRuns
             Use_Events(run,:)={[SplitConditions{run,1},ConfoundEvents{run,1}],[SplitConditionNames;ConfoundEvents{run,2}]};
        end
    else
        for run=1:numRuns
            [SplitConditions,SplitConditionNames,~] = ComputeSplitConditions(Events{run,1},ConditionNames,'SplitVar',rand(size(Events{run,1},1),1),'TrialNum',TrialNum{run,1},'RandomSplit',1);
            Use_Events(run,:)={[SplitConditions,ConfoundEvents{run,1}],[SplitConditionNames;ConfoundEvents{run,2}]};
        end  
    end
    [Condition_tVals,brain_mask,~,CondReg_CorrMat] = GetActivationPatterns(...
        Use_Events,TimeCourses,ConfoundTCs,...
        boldTCs,boldMasks,SplitConditionNames,...
        'parGLM',parGLM,...
        'Compute_CondReg_CorrMat',Compute_CondReg_CorrMat,...
        'ContrastNames',ContrastNames,...
        'Normalize',Normalize,...
        'ConfoundTCsBySubj',ConfoundTCsBySubj,...
        'ResampleSizes',ResampleSizes,...
        'AM2',AM2,...
        'AM2_Events',AM2_Events,...
        'TimecourseShift',TimecourseShift,...
        'ParcellationMask',ParcellationMask,...
        'ResampleSizesBoldTC',ResampleSizesBoldTC,...
        'NormWithinRun',NormWithinRun,...
        'NormAcrossRun',NormAcrossRun,...
        'Contrasts',Contrasts,...
        'BatchSize',BatchSize,...
        'ResampleToFit',ResampleToFit,...
        'UseStat',UseStat);  

    if isempty(gmMask)
        use_gmMask=brain_mask;
    else
        use_gmMask=gmMask;
    end
    brainSize=size(brain_mask);
    ConfoundRSMs=cell(1);
    ConfoundRSMLabels=cell(1);
    if Compute_CondReg_CorrMat==1
        ConfoundRSMs=cat(1,ConfoundRSMs,{CondReg_CorrMat});
        ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'CondReg_CorrMat'});
    end    
    if ~isempty(csfMask)
        ParcellationVector=csfMask(brain_mask~=0); 
        [ CSF_RSM ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
        ConfoundRSMs=cat(1,ConfoundRSMs,{CSF_RSM});
        ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'CSF_RSM'});       
    else
        CSF_RSM=[];
    end
    if ~isempty(wmMask)
        ParcellationVector=wmMask(brain_mask~=0); 
        [ WM_RSM ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
        ConfoundRSMs=cat(1,ConfoundRSMs,{WM_RSM});
        ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'WM_RSM'});        
    else
        WM_RSM=[];
    end    
    if ~isempty(ConfoundRSMs)
        if isempty(ConfoundRSMs{1,1})
            ConfoundRSMs(1,:)=[];
            ConfoundRSMLabels(1,:)=[];
        end
        numConfoundRSMs=size(ConfoundRSMs,1);
    end
    
    RSMs=cell(numParcels,1);
    
    if OutputConfoundRSMs==1
        All_ConfoundRSMs=cell2table(ConfoundRSMs,'VariableNames',ConfoundRSMLabels);
        if Compute_MeanAct_RSM==1
            All_ConfoundRSMs.MeanActRSMs=cell(numParcels,1);
        end
    else
        All_ConfoundRSMs=[];
    end
    RF=cell(numParcels,1);
    RC=cell(numParcels,1);
    IF=cell(numParcels,1);
    RFs=cell(numParcels,1);
    RCs=cell(numParcels,1);
    RFx=cell(numParcels,1);
    RCx=cell(numParcels,1);
    for i = 1:numParcels          
        UseMask=Parcels{i,1}.UseMask;
        UseLabels=Parcels{i,1}.UseLabels(:);
        if ~ismember(brainSize,size(UseMask),'rows')
            UseMask=imresize3(UseMask,brainSize,'nearest');
        end    
        gmVector=use_gmMask(brain_mask~=0);
        ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
        rsm_mask{i,rep}=single(brain_mask~=0).*UseMask.*use_gmMask;
        [ RSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
        if Compute_MeanAct_RSM==1
            [ MeanActRSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,'meanSim');
            if isempty(ConfoundRSMs)
                ConfoundRSMs{1,1}=MeanActRSMs;
            else
                ConfoundRSMs{numConfoundRSMs+1,1}=MeanActRSMs;
                ConfoundRSMLabels{numConfoundRSMs+1,1}='MeanActRSMs';
            end
            if OutputConfoundRSMs==1
                All_ConfoundRSMs.MeanActRSMs{i,1}=MeanActRSMs;
            end
        end
        if ~isempty(ConfoundRSMs)
            RSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs);
        end
        [RF{i,1},RC{i,1},IF{i,1},RFs{i,1},RFx{i,1},RCs{i,1},RCx{i,1}] = SplitHalfRCA(RSMs);
    end 
    Temp_RC{rep,1}=RC;
    Temp_RF{rep,1}=RF;
    Temp_IF{rep,1}=IF;
    Temp_RCs{rep,1}=RCs;
    Temp_RFs{rep,1}=RFs;
    Temp_RCx{rep,1}=RCx;
    Temp_RFx{rep,1}=RFx;    
end
rsm_mask=rsm_mask(:,1);
RC=cell(1,numParcels);
RF=cell(1,numParcels);
IF=cell(1,numParcels);
RCs=cell(1,numParcels);
RCx=cell(1,numParcels);
RFs=cell(1,numParcels);
RFx=cell(1,numParcels);
for i = 1:numParcels  
    for j = 1:NumReps
        All_RC{j,1}=Temp_RC{j,1}{i,1};
        All_RF{j,1}=Temp_RF{j,1}{i,1};
        All_IF{j,1}=Temp_IF{j,1}{i,1};
        All_RCs{j,1}=Temp_RCs{j,1}{i,1};
        All_RCx{j,1}=Temp_RCx{j,1}{i,1};
        All_RFs{j,1}=Temp_RFs{j,1}{i,1};
        All_RFx{j,1}=Temp_RFx{j,1}{i,1};
    end
    RC{1,i}=real(squeeze(cell2nDMAT(All_RC))');
    RF{1,i}=real(squeeze(cell2nDMAT(All_RF))');
    IF{1,i}=real(squeeze(cell2nDMAT(All_IF))');
    RCs{1,i}=real(squeeze(cell2nDMAT(All_RCs))');
    RCx{1,i}=real(squeeze(cell2nDMAT(All_RCx))');
    RFs{1,i}=real(squeeze(cell2nDMAT(All_RFs))');
    RFx{1,i}=real(squeeze(cell2nDMAT(All_RFx))');
end 

end

function CleanRSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs)
    numConfoundRSMs=size(ConfoundRSMs,1);
    RSMSize=size(RSMs,1);
    RSMs = mat2uppertriuvectormat(RSMs);
    numRSMs=size(RSMs,2);
    indMat=ones(numRSMs,numConfoundRSMs);
    residRSMs=RSMs*0;
    RegConstant=ones(size(RSMs,1),1);
    for i = 1:numConfoundRSMs
        ConfoundRSMs{i,1}=mat2uppertriuvectormat(ConfoundRSMs{i,1});
        if size(ConfoundRSMs{i,1},2)==numRSMs
            indMat(:,i)=[1:numRSMs]';
        end
    end
    parfor i = 1:numRSMs
        tempRSMs=RSMs(:,i);
        tempConfoundRSMs=RegConstant;
        for j = 1:numConfoundRSMs
            tempConfoundRSMs=[tempConfoundRSMs,ConfoundRSMs{j,1}(:,indMat(i,j))];
        end
        [residRSMs(:,i)] = FastOLSRegress_Resids(tempRSMs,tempConfoundRSMs);
    end
    [ CleanRSMs ] = vertRSM2SymRSM( residRSMs,RSMSize );    
end
