function [All_RSMs,rsm_mask,SplitData,DesignMatrix,AfniInfo,Condition_tVals1,Condition_tVals2,Condition_bVals1,Condition_bVals2,slRSMs,SearchlightResults] = GetParcellation2SplitVarRCA(Events,ConfoundEvents,TimeCourses,ConfoundTCs,boldTCs,boldMasks,ConditionNames,Parcels,TrialNum,SplitVar,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ActivationPatterns] = VariableSetter('ActivationPatterns',[],varargin);
[numGLM] = VariableSetter('numGLM',2,varargin);
[UseAfni] = VariableSetter('UseAfni',0,varargin);
[TrialOnsetTimes] = VariableSetter('TrialOnsetTimes',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
[UseAfniConfounds] = VariableSetter('UseAfniConfounds',1,varargin);
Use3DReml = VariableSetter('Use3DReml',0,varargin);
[SkipRSMs] = VariableSetter('SkipRSMs',0,varargin);
[TR] = VariableSetter('TR',[],varargin);
[SearchlightInfo] = VariableSetter('SearchlightInfo',[],varargin);
[VoxNorm] = VariableSetter('VoxNorm',[],varargin);
%%% Optional input for GetParcellationSplitRCA
[SplitAcrossRun] = VariableSetter('SplitAcrossRun',0,varargin);
%%% Optional input for RSMParcellationRSM
[csfMask] = VariableSetter('csfMask',[],varargin);
[wmMask] = VariableSetter('wmMask',[],varargin);
[invgmMask] = VariableSetter('invgmMask',[],varargin);
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
Condition_bVals1=[];
Condition_bVals2=[];
All_RSMs=[];
SplitData=[];
DesignMatrix=[];
AfniInfo=[];
Condition_tVals1=[];
Condition_tVals2=[];
numRuns=size(Events,1);
numParcels=size(Parcels,1);
SearchlightResults=SearchlightInfo;
slRSMs=[];
if isempty(ActivationPatterns)
    if SplitAcrossRun==1
        [EventsMat,CatInd] = Cell2CatMat(Events(:,1));
        [SplitVar] = Cell2CatMat(SplitVar);
        TrialMat=[];
        for run=1:numRuns
            if run > 1
                MaxVal=max(TrialMat(:));
                TrialMat=[TrialMat;TrialNum{run,1}+MaxVal];
            else
                TrialMat=TrialNum{run,1};
            end
        end
        SplitData=zeros(1,2);
    else
        SplitData=zeros(numRuns,2);
        EventsMat=[];
        TrialMat=[];
        CatInd=[];
    end   
    if numGLM == 2
        Use_Events1=cell(numRuns,2);
        Use_Events2=cell(numRuns,2);
        if SplitAcrossRun==1
            [SplitConditions,~,~,tempSplitData] = ComputeSplitConditions(EventsMat,ConditionNames,'SplitVar',SplitVar,'TrialNum',TrialMat);
            SplitData=tempSplitData.groupMeans;
            SplitConditions1=SplitConditions(:,1:(size(SplitConditions,2)/2));
            SplitConditions2=SplitConditions(:,(size(SplitConditions,2)/2)+1:end);
            Split1_Events=single(nansum(SplitConditions1,2)>0); 
            Split2_Events=single(nansum(SplitConditions2,2)>0);     
            [SplitConditions1] = CatMat2Cell([SplitConditions1,Split2_Events],CatInd); 
            [SplitConditions2] = CatMat2Cell([SplitConditions2,Split1_Events],CatInd);
            for run=1:numRuns
                Use_Events1(run,:)={[SplitConditions1{run,1},ConfoundEvents{run,1}],[ConditionNames;{'Split2_Events'};ConfoundEvents{run,2}]};
                Use_Events2(run,:)={[SplitConditions2{run,1},ConfoundEvents{run,1}],[ConditionNames;{'Split1_Events'};ConfoundEvents{run,2}]};
            end
        else
            for run=1:numRuns
                [SplitConditions,~,~,tempSplitData] = ComputeSplitConditions(Events{run,1},ConditionNames,'SplitVar',SplitVar{run,1},'TrialNum',TrialNum{run,1});
                SplitData(run,:)=tempSplitData.groupMeans;
                SplitConditions1=SplitConditions(:,1:(size(SplitConditions,2)/2));
                SplitConditions2=SplitConditions(:,(size(SplitConditions,2)/2)+1:end); 
                Split1_Events=single(nansum(SplitConditions1,2)>0); 
                Split2_Events=single(nansum(SplitConditions2,2)>0);                       
                Use_Events1(run,:)={[SplitConditions1,Split2_Events,ConfoundEvents{run,1}],[ConditionNames;{'Split2_Events'};ConfoundEvents{run,2}]};
                Use_Events2(run,:)={[SplitConditions2,Split1_Events,ConfoundEvents{run,1}],[ConditionNames;{'Split1_Events'};ConfoundEvents{run,2}]};
            end  
            SplitData=nanmean(SplitData,1);
        end
        [Condition_tVals1,~,~,~,DesignMatrix1,AfniInfo1,Condition_bVals1] = GetActivationPatterns(...
            Use_Events1,TimeCourses,ConfoundTCs,...
            boldTCs,boldMasks,ConditionNames,...
            'parGLM',parGLM,...
            'UseAfni',UseAfni,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'UseAfniConfounds',UseAfniConfounds,...
            'Use3DReml',Use3DReml,...
            'AFNISkipDelete',1,...
            'SkipAfniInput',0,...
            'TR',TR,...
            'TrialNum',TrialNum,...        
            'Compute_CondReg_CorrMat',0,...
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
        [Condition_tVals2,brain_mask,~,~,DesignMatrix2,AfniInfo2,Condition_bVals2] = GetActivationPatterns(...
            Use_Events2,TimeCourses,ConfoundTCs,...
            boldTCs,boldMasks,ConditionNames,...
            'parGLM',parGLM,...
            'UseAfni',UseAfni,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'UseAfniConfounds',UseAfniConfounds,...
            'Use3DReml',Use3DReml,...
            'AFNISkipDelete',0,...
            'SkipAfniInput',1,...        
            'TR',TR,...
            'TrialNum',TrialNum,...        
            'Compute_CondReg_CorrMat',0,...
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

            AfniInfo={AfniInfo1,AfniInfo2};
            DesignMatrix={DesignMatrix1,DesignMatrix2};
                
    else
        Use_Events=cell(numRuns,2);
        if SplitAcrossRun==1
            [SplitConditions,SplitConditionNames,~,tempSplitData] = ComputeSplitConditions(EventsMat,ConditionNames,'SplitVar',SplitVar,'TrialNum',TrialMat);
            SplitData=tempSplitData.groupMeans;
            [SplitConditions] = CatMat2Cell(SplitConditions,CatInd);        
            for run=1:numRuns
                 Use_Events(run,:)={[SplitConditions{run,1},ConfoundEvents{run,1}],[SplitConditionNames;ConfoundEvents{run,2}]};
            end
        else
            for run=1:numRuns
                [SplitConditions,SplitConditionNames,~,tempSplitData] = ComputeSplitConditions(Events{run,1},ConditionNames,'SplitVar',SplitVar{run,1},'TrialNum',TrialNum{run,1});
                SplitData(run,:)=tempSplitData.groupMeans;
                if ~isempty(ConfoundEvents{1,1})
                    Use_Events(run,:)={[SplitConditions,ConfoundEvents{run,1}],[SplitConditionNames;ConfoundEvents{run,2}]};
                end
            end  
            SplitData=nanmean(SplitData,1);
        end
        [Condition_tVals,brain_mask,~,~,DesignMatrix,AfniInfo,Condition_bVals] = GetActivationPatterns(...
            Use_Events,TimeCourses,ConfoundTCs,...
            boldTCs,boldMasks,SplitConditionNames,...
            'parGLM',parGLM,...
            'UseAfni',UseAfni,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'UseAfniConfounds',UseAfniConfounds,...
            'Use3DReml',Use3DReml,...
            'AFNISkipDelete',0,...
            'SkipAfniInput',0,...
            'TR',TR,...
            'TrialNum',TrialNum,...        
            'Compute_CondReg_CorrMat',0,...
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
        Condition_tVals1=Condition_tVals(:,1:size(Condition_tVals,2)/2);
        Condition_tVals2=Condition_tVals(:,(1+size(Condition_tVals,2)/2):end);
        Condition_bVals1=Condition_bVals(:,1:size(Condition_bVals,2)/2);
        Condition_bVals2=Condition_bVals(:,(1+size(Condition_bVals,2)/2):end); 
        DesignMatrix1=DesignMatrix;
        DesignMatrix2=DesignMatrix;
        DesignMatrix1(:,SplitConditionNames((1+size(Condition_bVals,2)/2):end,1))=[];
        DesignMatrix1.Properties.VariableNames(ismember(DesignMatrix1.Properties.VariableNames,SplitConditionNames))=ConditionNames;
        DesignMatrix2(:,SplitConditionNames(1:size(Condition_bVals,2)/2,1))=[];
        DesignMatrix2.Properties.VariableNames(ismember(DesignMatrix2.Properties.VariableNames,SplitConditionNames))=ConditionNames;
    end    
else
    Condition_tVals1=ActivationPatterns.Split1_ActVals;
    Condition_tVals2=ActivationPatterns.Split2_ActVals;
    brain_mask=ActivationPatterns.brain_mask;
    AfniInfo=ActivationPatterns.AfniInfo;
    DesignMatrix=ActivationPatterns.DesignMatrix;
    ConditionNames=ActivationPatterns.AnalysisParameters.ConditionNames;
    SplitData=ActivationPatterns.SplitData;
    if iscell(DesignMatrix)
        DesignMatrix1=DesignMatrix{1,1};
        DesignMatrix2=DesignMatrix{1,2};
    elseif istable(DesignMatrix)
        SplitConditionNames1=join([ConditionNames,repmat({'_1'},[length(ConditionNames),1])],'');
        SplitConditionNames2=join([ConditionNames,repmat({'_2'},[length(ConditionNames),1])],'');
        DesignMatrix1=DesignMatrix;
        DesignMatrix2=DesignMatrix;
        DesignMatrix1(:,SplitConditionNames2)=[];
        DesignMatrix1.Properties.VariableNames(ismember(DesignMatrix1.Properties.VariableNames,SplitConditionNames1))=ConditionNames;
        DesignMatrix2(:,SplitConditionNames1)=[];
        DesignMatrix2.Properties.VariableNames(ismember(DesignMatrix2.Properties.VariableNames,SplitConditionNames2))=ConditionNames;   
    end
    
end
   
if SkipRSMs==0
    %All_RSMs=zeros((size(Condition_tVals1,1)^2-1)/2,2,'single');   
    if strcmpi(UseStat,'B') && isempty(ActivationPatterns)
        Condition_tVals2=Condition_bVals2;
        Condition_tVals1=Condition_bVals1;    
        Condition_bVals1=[];
        Condition_bVals2=[];
    end

    Condition_tVals=Condition_tVals1;
    if strcmpi(VoxNorm,'mean')
        Condition_tVals=Condition_tVals-repmat(nanmean(Condition_tVals,2),[1,size(Condition_tVals,2)]);
    elseif strcmpi(VoxNorm,'Z')
        Condition_tVals=nan_zscore(Condition_tVals')';
    end
        
    if isempty(gmMask)
        use_gmMask=brain_mask;
    else
        use_gmMask=gmMask;
    end
    All_RSMs=cell(numParcels,1);
    brainSize=size(brain_mask);
    ConfoundRSMs=cell(1);
    ConfoundRSMLabels=cell(1);
    if Compute_CondReg_CorrMat==1
        CondReg_CorrMat=corrcoef(table2array(DesignMatrix1(:,ConditionNames)));
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
    
    if OutputConfoundRSMs==1
        All_ConfoundRSMs=cell2table(ConfoundRSMs,'VariableNames',ConfoundRSMLabels);
        if Compute_MeanAct_RSM==1
            All_ConfoundRSMs.MeanActRSMs=cell(numParcels,1);
        end
    else
        All_ConfoundRSMs=[];
    end
     rsm_mask=cell(1);
    for i = 1:numParcels    
        if strcmpi(Parcels{i,1},'Searchlight')
            rsm_mask=brain_mask;
             [ RSMs,MeanActRSMs,rsm_mask,SearchlightResults.SLinds,SearchlightResults.numVoxels,SearchlightResults.MaskCoords] = SearchlightRSM( Condition_tVals,brain_mask,SearchlightInfo.slRadius,SearchlightInfo.slVoxelThreshold,SearchlightInfo.slShape,1000,...
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
                RSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs);
            end     

        else           
            UseMask=Parcels{i,1}.UseMask;
            UseLabels=Parcels{i,1}.UseLabels(:);
            if ~ismember(brainSize,size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize,'nearest');
            end    
            gmVector=use_gmMask(brain_mask~=0);
            ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
            [ RSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
            RSMs = mat2uppertriuvectormat(RSMs);
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
        end
        All_RSMs{i,1}=RSMs;
    end 
    
    Condition_tVals=Condition_tVals2;
    if strcmpi(VoxNorm,'mean')
        Condition_tVals=Condition_tVals-repmat(nanmean(Condition_tVals,2),[1,size(Condition_tVals,2)]);
    elseif strcmpi(VoxNorm,'Z')
        Condition_tVals=nan_zscore(Condition_tVals')';
    end    
    if isempty(gmMask)
        use_gmMask=brain_mask;
    else
        use_gmMask=gmMask;
    end
    brainSize=size(brain_mask);
    ConfoundRSMs=cell(1);
    ConfoundRSMLabels=cell(1);
    if Compute_CondReg_CorrMat==1
        CondReg_CorrMat=corrcoef(table2array(DesignMatrix2(:,ConditionNames)));
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
             [ RSMs,MeanActRSMs,rsm_mask,SearchlightResults.SLinds,SearchlightResults.numVoxels,SearchlightResults.MaskCoords] = SearchlightRSM( Condition_tVals,brain_mask,SearchlightInfo.slRadius,SearchlightInfo.slVoxelThreshold,SearchlightInfo.slShape,100,...
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
                RSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs);
            end     

        else            
            UseMask=Parcels{i,1}.UseMask;
            UseLabels=Parcels{i,1}.UseLabels(:);
            if ~ismember(brainSize,size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize,'nearest');
            end    
            gmVector=use_gmMask(brain_mask~=0);
            ParcellationVector=UseMask(brain_mask~=0).*gmVector;            
            rsm_mask{i,1}=single(brain_mask~=0).*UseMask.*use_gmMask;
            [ RSMs ] = ParcellationRSMs( Condition_tVals,ParcellationVector,SimType);
            RSMs = mat2uppertriuvectormat(RSMs);
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
        end
        All_RSMs{i,2}=RSMs;
    end 
else
    rsm_mask=brain_mask;
end

end

function CleanRSMs = ComputeCleanRSMs(RSMs,ConfoundRSMs)
    numConfoundRSMs=size(ConfoundRSMs,1);
    numRSMs=size(RSMs,2);
    indMat=ones(numRSMs,numConfoundRSMs);
    residRSMs=RSMs*0;
    %RegConstant=ones(size(RSMs,1),1);
    for i = 1:numConfoundRSMs
        ConfoundRSMs{i,1}=mat2uppertriuvectormat(ConfoundRSMs{i,1});
        if size(ConfoundRSMs{i,1},2)==numRSMs
            indMat(:,i)=[1:numRSMs]';
        end
    end
    parfor i = 1:numRSMs
        tempRSMs=RSMs(:,i);
        tempConfoundRSMs=[];
       %tempConfoundRSMs=RegConstant;
        for j = 1:numConfoundRSMs
            tempConfoundRSMs=[tempConfoundRSMs,ConfoundRSMs{j,1}(:,indMat(i,j))];
        end
        tempstats=regstats(tempRSMs,tempConfoundRSMs,'linear','standres');
        residRSMs(:,i)=tempstats.standres;
        %[residRSMs(:,i)] = FastOLSRegress_Resids(tempRSMs,tempConfoundRSMs);
    end
    %[ CleanRSMs ] = vertRSM2SymRSM( residRSMs,RSMSize );    
    CleanRSMs=residRSMs;
end
