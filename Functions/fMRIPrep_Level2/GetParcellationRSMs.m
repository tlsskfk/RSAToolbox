function [RSMs,rsm_mask,All_ConfoundRSMs,RegressorNames,AfniInfo,slRSMs,SearchlightResults] = GetParcellationRSMs(Events,TimeCourses,ConfoundTCs,boldTCs,boldMasks,ConditionNames,Parcels,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ActivationPatterns] = VariableSetter('ActivationPatterns',[],varargin);
[UseAfni] = VariableSetter('UseAfni',0,varargin);
Use3DReml = VariableSetter('Use3DReml',0,varargin);
[TrialOnsetTimes] = VariableSetter('TrialOnsetTimes',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
[TR] = VariableSetter('TR',[],varargin);
[TrialNum] = VariableSetter('TrialNum',[],varargin);
%%% Optional input for RSMParcellationRSM
[csfMask] = VariableSetter('csfMask',[],varargin);
[wmMask] = VariableSetter('wmMask',[],varargin);
[invgmMask] = VariableSetter('invgmMask',[],varargin);
[Compute_MeanAct_RSM] = VariableSetter('Compute_MeanAct_RSM',0,varargin);
[gmMask] = VariableSetter('gmMask',[],varargin);
[SimType] = VariableSetter('SimType','corrcoef',varargin);
[OutputConfoundRSMs] = VariableSetter('OutputConfoundRSMs',0,varargin);
[SearchlightInfo] = VariableSetter('SearchlightInfo',[],varargin);
[VoxNorm] = VariableSetter('VoxNorm',[],varargin);
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

SearchlightResults=SearchlightInfo;
slRSMs=[];
if isempty(ActivationPatterns)
    [Condition_tVals,brain_mask,RegressorNames,CondReg_CorrMat,~,AfniInfo] = GetActivationPatterns(...
        Events,TimeCourses,ConfoundTCs,...
        boldTCs,boldMasks,ConditionNames,...
        'parGLM',parGLM,...
        'UseAfni',UseAfni,...
        'Use3DReml',Use3DReml,...
        'TrialOnsetTimes',TrialOnsetTimes,...
        'AfniWorkDir',AfniWorkDir,...
        'TR',TR,...
        'TrialNum',TrialNum,...    
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
else
    Condition_tVals=ActivationPatterns.ActVals;
    brain_mask=ActivationPatterns.brain_mask;
    try
        AfniInfo=ActivationPatterns.AfniInfo;
        DesignMatrix=ActivationPatterns.DesignMatrix;
        RegressorNames=DesignMatrix.Properties.VariableNames;
        if Compute_CondReg_CorrMat==1
            ConditionNames=ActivationPatterns.AnalysisParameters.ConditionNames;
            CondReg_CorrMat=corrcoef(table2array(DesignMatrix(:,ConditionNames)));
        else
            CondReg_CorrMat=[];
        end
    catch
        CondReg_CorrMat=[];
        RegressorNames=[];
        AfniInfo=[];
    end
end
if isempty(gmMask)
    gmMask=brain_mask;
end

brainSize=size(brain_mask);
ConfoundRSMs=cell(1);
ConfoundRSMLabels=cell(1);
if strcmpi(VoxNorm,'mean')
    Condition_tVals=Condition_tVals-repmat(nanmean(Condition_tVals,2),[1,size(Condition_tVals,2)]);
elseif strcmpi(VoxNorm,'Z')
    Condition_tVals=nan_zscore(Condition_tVals')';
end
    
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
if ~isempty(invgmMask)
    ParcellationVector=invgmMask(brain_mask~=0); 
    [ GM_RSM ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
    ConfoundRSMs=cat(1,ConfoundRSMs,{GM_RSM});
    ConfoundRSMLabels=cat(1,ConfoundRSMLabels,{'GM_RSM'});        
else
    GM_RSM=[];
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
    if strcmpi(Parcels{i,1},'Searchlight')
         [ RSMs{i,1},MeanActRSMs,rsm_mask{i,1},SearchlightResults.SLinds,SearchlightResults.numVoxels,SearchlightResults.MaskCoords] = SearchlightRSM( Condition_tVals,brain_mask,SearchlightInfo.slRadius,SearchlightInfo.slVoxelThreshold,SearchlightInfo.slShape,100,...
             'SLinds',SearchlightInfo.SLinds,...
             'numVoxels',SearchlightInfo.SLnumVoxels,...
             'MaskCoords',SearchlightInfo.SLcoords,...
             'ComputeMeanRSM',Compute_MeanAct_RSM);
        if Compute_MeanAct_RSM==1
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
        
    else    
        UseMask=Parcels{i,1}.UseMask;
        UseLabels=Parcels{i,1}.UseLabels(:);
        if ~ismember(brainSize,size(UseMask),'rows')
            UseMask=imresize3(UseMask,brainSize,'nearest');
        end    
        gmVector=gmMask(brain_mask~=0);
        ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
        rsm_mask{i,1}=single(brain_mask~=0).*UseMask.*gmMask;
        [ RSMs{i,1} ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType,length(UseLabels));
        RSMs{i,1} = mat2uppertriuvectormat(RSMs{i,1});
        if Compute_MeanAct_RSM==1
            [ MeanActRSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,'meanSim',length(UseLabels));
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
    end
    if ~isempty(ConfoundRSMs)
        RSMs{i,1} = ComputeCleanRSMs(RSMs{i,1},ConfoundRSMs);
    end               
end 

end

function CleanRSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs)
    numConfoundRSMs=size(ConfoundRSMs,1);
    %RSMSize=size(RSMs,1);
    %RSMs = mat2uppertriuvectormat(RSMs);
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
    for i = 1:numRSMs
        tempRSMs=RSMs(:,i);
        tempConfoundRSMs=[];
        for j = 1:numConfoundRSMs
            tempConfoundRSMs=[tempConfoundRSMs,ConfoundRSMs{j,1}(:,indMat(i,j))];
        end
        try
            tempstats=regstats(tempRSMs,tempConfoundRSMs,'linear','standres');
            residRSMs(:,i)=tempstats.standres;
        catch
            residRSMs(:,i)=nan(length(RegConstant),1);
        end
    end
    %[ CleanRSMs ] = vertRSM2SymRSM( residRSMs,RSMSize );    
    CleanRSMs=residRSMs;
end
