function [RSMs,rsm_mask,All_ConfoundRSMs,RegressorNames] = GetParcellationSplitRCA(Events,TimeCourses,ConfoundTCs,boldTCs,boldMasks,ConditionNames,Parcels,SplitParams,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

parfor 

[Condition_tVals,brain_mask,RegressorNames,CondReg_CorrMat] = GetActivationPatterns(...
    Events,TimeCourses,ConfoundTCs,...
    boldTCs,boldMasks,ConditionNames,...
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
    gmMask=brain_mask;
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
    ParcellationVector=Use_wmMask(brain_mask~=0); 
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
numParcels=size(Parcels,1);
RSMs=cell(numParcels,1);
rsm_mask=cell(numParcels,1);
if OutputConfoundRSMs==1
    All_ConfoundRSMs=cell2table(ConfoundRSMs,'VariableNames',ConfoundRSMLabels);
    if Compute_MeanAct_RSM==1
        All_ConfoundRSMs.MeanActRSMs=cell(numParcels,1);
    end
else
    All_ConfoundRSMs=[];
end
        
for i = 1:numParcels          
    UseMask=Parcels{i,1}.UseMask;
    UseLabels=Parcels{i,1}.UseLabels(:);
    if ~ismember(brainSize,size(UseMask),'rows')
        UseMask=imresize3(UseMask,brainSize,'nearest');
    end    
    gmVector=gmMask(brain_mask~=0);
    ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
    rsm_mask{i,1}=single(brain_mask~=0).*UseMask.*gmMask;
    [ RSMs{i,1} ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
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
        RSMs{i,1} = ComputeCleanRSMs(RSMs{i,1},ConfoundRSMs);
    end               
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
