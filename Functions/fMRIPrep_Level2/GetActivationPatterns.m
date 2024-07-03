function [Condition_tVals,brain_mask,AllRegressorNames,CondReg_CorrMat,DesignMatrix,AfniInfo,Condition_bVals] = GetActivationPatterns(Events,TimeCourses,ConfoundTCs,boldTCs,boldMasks,ConditionNames,varargin)
%Written by David Rothlein
%%% Optional input for GetActivationPatterns
[UseAfni] = VariableSetter('UseAfni',0,varargin);
Use3DReml = VariableSetter('Use3DReml',0,varargin);
[TrialOnsetTimes] = VariableSetter('TrialOnsetTimes',[],varargin);
[AfniWorkDir] = VariableSetter('AfniWorkDir',[],varargin);
AFNISkipDelete = VariableSetter('AFNISkipDelete',0,varargin);
[TR] = VariableSetter('TR',[],varargin);
[TrialNum] = VariableSetter('TrialNum',[],varargin);
[SkipAfniInput] = VariableSetter('SkipAfniInput',0,varargin);
[parGLM] = VariableSetter('parGLM',0,varargin);
[Compute_CondReg_CorrMat] = VariableSetter('Compute_CondReg_CorrMat',0,varargin);
[ContrastNames] = VariableSetter('ContrastNames',[],varargin);
[UseAfniConfounds] = VariableSetter('UseAfniConfounds',1,varargin);
%%% Optional input for GLM_RegressorPrep
[Normalize] = VariableSetter('Normalize','zscore',varargin);
[ConfoundTCsBySubj] = VariableSetter('ConfoundTCsBySubj',{[]},varargin);
[ResampleSizes] = VariableSetter('ResampleSizes',[],varargin);
[AM2]=VariableSetter('AM2',1,varargin);
[AM2_Events]=VariableSetter('AM2_Events',{[]},varargin);
TimecourseShift = VariableSetter('TimecourseShift',0,varargin);
[VoxNorm] = VariableSetter('VoxNorm',[],varargin);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AfniInfo=[];
if UseAfni==0
    numRuns=size(boldTCs,1);
    if ~isempty(ConfoundTCs) || ~isempty(ConfoundTCsBySubj)
        for run = 1:numRuns
            if ~isempty(ConfoundTCsBySubj) 
                if size(ConfoundTCsBySubj{run,1},1)>ResampleSizes{run,1}
                    ConfoundTCsBySubj{run,1}=ConfoundTCsBySubj{run,1}(1:ResampleSizes{run,1},:);
                end
            end
            if ~isempty(ConfoundTCs) 
                if size(ConfoundTCs{run,1},1)>ResampleSizes{run,1}
                    ConfoundTCs{run,1}=ConfoundTCs{run,1}(1:ResampleSizes{run,1},:);
                end
            end   
        end
    end
    [Regressors,RegressorNames,~,EmptyRegressors] = GLM_RegressorPrep(Events,...
        TimeCourses,...
        ConfoundTCs,...
        'AM2',AM2,...
        'AM2_Events',AM2_Events,...
        'TimecourseShift',TimecourseShift,...
        'Normalize',Normalize,...
        'ConfoundTCsBySubj',ConfoundTCsBySubj,...
        'ResampleSizes',ResampleSizes);  

    AllRegressorNames=[RegressorNames;ContrastNames];

    [BoldTC,brain_mask] = GLM_boldTCPrep(boldTCs,boldMasks,...
        'ParcellationMask',ParcellationMask,...
        'ResampleSizesBoldTC',ResampleSizesBoldTC,...
        'NormWithinRun',NormWithinRun,...
        'NormAcrossRun',NormAcrossRun);

    UseInd=ismember(AllRegressorNames,ConditionNames);
    RegressorNames=RegressorNames(EmptyRegressors==0);
    if parGLM==1
        [Condition_tVals,DesignMatrix,Condition_bVals] = FastGLM_parfor(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);    
    else    
        [Condition_tVals,DesignMatrix,Condition_bVals] = FastGLM(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);
    end
    DesignMatrix=array2table(DesignMatrix,'VariableNames',RegressorNames);
    CondReg_CorrMat=[];
    if Compute_CondReg_CorrMat==1
        CondReg_CorrMat=corrcoef(Regressors(:,UseInd==1));
    end
    Condition_tVals=Condition_tVals(:,UseInd==1);
    Condition_bVals=Condition_bVals(:,UseInd==1);     
elseif UseAfni==1
    numRuns=size(boldTCs,1);
    jobs = feature('numcores');
    if ispc
        if strcmp(AfniWorkDir(1,2),':')
            DirLetter=lower(AfniWorkDir(1,1));
            AltAfniWorkDir=strrep(AfniWorkDir,AfniWorkDir(1,1:2),['/mnt/',DirLetter]);
        end
    else
        AltAfniWorkDir=AfniWorkDir;
    end    

    % Create timecourse files
    if any(~cellfun(@isempty,ConfoundTCs(:,1)))    
        for run = 1:numRuns
            if isempty(ConfoundTCs{run,1})
                useInd=find(~cellfun(@isempty,ConfoundTCs(:,1)));
                ConfoundTCs{run,1}=zeros(ResampleSizes{run,1},size(ConfoundTCs{useInd(1,1),1},2));
                ConfoundTCs{run,2}=ConfoundTCs{useInd(1,1),2};
            end
        end
    end
    if iscell(ConfoundTCs) && isempty(ConfoundTCs{1,1})
        ConfoundTCs=[];
    end
    if iscell(ConfoundTCsBySubj) && isempty(ConfoundTCsBySubj{1,1})
        ConfoundTCsBySubj=[];
    end    
    if ~isempty(ConfoundTCs) || ~isempty(ConfoundTCsBySubj)
        for run = 1:numRuns
            if ~isempty(ConfoundTCsBySubj) 
                if size(ConfoundTCsBySubj{run,1},1)>ResampleSizes{run,1}
                    ConfoundTCsBySubj{run,1}=ConfoundTCsBySubj{run,1}(1:ResampleSizes{run,1},:);
                end
            end
            if ~isempty(ConfoundTCs) 
                if size(ConfoundTCs{run,1},1)>ResampleSizes{run,1}
                    ConfoundTCs{run,1}=ConfoundTCs{run,1}(1:ResampleSizes{run,1},:);
                end
            end   
        end
        [tempRegressors,tempRegressorNames,~,~,tempRegressorNamesRaw] = GLM_RegressorPrep(Events,...
            TimeCourses,...
            ConfoundTCs,...
            'AM2',AM2,...
            'AM2_Events',AM2_Events,...
            'TimecourseShift',TimecourseShift,...
            'Normalize',Normalize,...
            'ConfoundTCsBySubj',ConfoundTCsBySubj,...
            'ResampleSizes',ResampleSizes);          
        AllConfoundTCNames=[];
        tempCell=cell(numRuns,1);

        for run=1:numRuns
            if ~isempty(ConfoundTCsBySubj)
                AllConfoundTCNames=unique([AllConfoundTCNames;ConfoundTCsBySubj{run,2}]);
                tempCell{run,1}=ConfoundTCsBySubj{run,1}(:,1);
            end
            if ~isempty(ConfoundTCs)
                AllConfoundTCNames=unique([AllConfoundTCNames;ConfoundTCs{run,2}]);
                tempCell{run,1}=ConfoundTCs{run,1}(:,1);
            end          
        end
        [~,CatInd] = Cell2CatMat(tempCell);
        ConfoundTCs=[CatMat2Cell(tempRegressors(:,ismember(tempRegressorNamesRaw,AllConfoundTCNames)),CatInd),repmat({tempRegressorNames(ismember(tempRegressorNamesRaw,AllConfoundTCNames),:)},[numRuns,1])];
        numConfoundTCs=size(ConfoundTCs{1,1},2);
        BaseLine_Info=cell(numConfoundTCs,2);
        for ConfoundTCNum = 1:numConfoundTCs
            tempConfoundTC = cell(numRuns,1);
            for run=1:numRuns
                try
                    tempConfoundTC{run,1} = ConfoundTCs{run,1}(:,ConfoundTCNum);
                catch
                    tempConfoundTC{run,1} = zeros(size(ConfoundTCs{run,1},1),1);
                end
            end
            BaseLine_Info{ConfoundTCNum,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
%             AFNILabels{AFNILabelCount,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
%             AFNILabelCount=AFNILabelCount+1;
            [BaseLine_Info{ConfoundTCNum,2}] = AFNIonsetMaker_stimfile(tempConfoundTC,AfniWorkDir,['ConfoundTC',num2str(ConfoundTCNum)],AltAfniWorkDir);
        end
    else
        BaseLine_Info=[];
    end   
    AFNILabels=cell(1);
    AFNILabelCount=1;    
    if ~isempty(Events)
        numEvents=size(Events{1,1},2); 
        Event_Info=cell(numEvents,2);
        for eventNum = 1:numEvents
            tempEvents=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempEvents{run,1}=Events{run,1}(tInd,eventNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
            end
            Event_Info{eventNum,1}=Events{1,2}{eventNum,1};
            AFNILabels{AFNILabelCount,1}=Events{1,2}{eventNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [Event_Info{eventNum,2}] = AFNIonsetMaker_event(tempEvents,tempOnsets,AfniWorkDir,['Event',num2str(eventNum)],AltAfniWorkDir);
        end
    else
        Event_Info=[];
    end
%     if ~isempty(AM2_Events)
%         numEvents=size(AM2_Events{1,1},2); 
%         AM2_Events_Info=cell(numEvents,2);
%         for eventNum = 1:numEvents
%             tempEvents=cell(numRuns,1);
%             tempOnsets=cell(numRuns,1);
%             for run=1:numRuns
%                 [~,tInd]=unique(TrialNum{run,1});
%                 tempEvents{run,1}=AM2_Events{run,1}(tInd,eventNum);
%                 tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
%             end
%             AM2_Events_Info{eventNum,1}=AM2_Events{1,2}{eventNum,1};
%             AFNILabels{AFNILabelCount,1}=AM2_Events{1,2}{eventNum,1};
%             AFNILabelCount=AFNILabelCount+1;
%             [AM2_Events_Info{eventNum,2}] = AFNIonsetMaker_event(tempEvents,tempOnsets,AfniWorkDir,['AM2_Event',num2str(eventNum)],AltAfniWorkDir);
%         end
%         Event_Info=[Event_Info;AM2_Events_Info];
%     else
%         AM2_Events_Info=[];
%     end  
    
    if ~isempty(TimeCourses{1,1})
        numTCs=size(TimeCourses{1,1},2);
        numAM2Events=size(AM2_Events{1,1},2);
        TC_Info=cell(numTCs,2);
        for tcNum = 1:numTCs
            tempTimeCourses=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            tempAM2_Events=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,tcNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
                if numAM2Events==0
                    tempAM2_Events{run,1}=ones(length(tempTimeCourses{run,1}),1);
                elseif numAM2Events==1
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
                else
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,tcNum);
                end
            end
            AFNILabels{AFNILabelCount,1}=['AM2_Event_',TimeCourses{1,2}{tcNum,1}];
            AFNILabelCount=AFNILabelCount+1;
            
            TC_Info{tcNum,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [TC_Info{tcNum,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['TCs',num2str(tcNum)],AltAfniWorkDir);
        end
    else
        TC_Info=[];
    end
    AllRegressorNames=AFNILabels;
    % Convert boldTC + mask to brik/head
    brain_mask=single(sum(single(cell2nDMAT(boldMasks)~=0),4)==size(boldMasks,1));
    inputNames=[];
    for run = 1:numRuns
        BrainSaveName=[AfniWorkDir,'Run',num2str(run)];
        AltBrainSaveName=[AltAfniWorkDir,'Run',num2str(run)];
        if SkipAfniInput==0
            [brainMap] = brainMask2brainMap(boldTCs{run,1}',boldMasks{run,1});
            brainMap=brainMap.*repmat(brain_mask,[1,1,1,size(brainMap,4)]);
            SaveBrik_3mmMNI(brainMap,[],BrainSaveName);
        end
        inputNames=[inputNames,' ',AltBrainSaveName,'+tlrc.HEAD'];
    end
% Generate AFNI Deconvolve script and run
    if Use3DReml==0
        [ runCode ] = AFNI_3dDeconvolve(inputNames,Event_Info,[AltAfniWorkDir,'glmMap'],...
            'polort','A',...
            'TR',TR,...
            'AM2_stimInfo',AM2TC_Info,...
            'BaseLine_stimInfo',BaseLine_Info,...
            'jobs',jobs,...
            'overwrite',1,...
            'HDR','GAM');
        AfniInfo.runCode=runCode;
        afniCodeFileName=[AfniWorkDir,'glmCode.tcsh'];
        ALTafniCodeFileName=[AltAfniWorkDir,'glmCode.tcsh'];
        fileID = fopen(afniCodeFileName,'w');
        fprintf(fileID,'%s',runCode);
        fclose(fileID);
        if ispc
            [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);
        else
            [~,AfniInfo.afniOut]=system(['tcsh ',ALTafniCodeFileName]);
        end

    % Load AFNI brik/head output into matlab and extract condition tVals/betas
        [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);

        UseDesignMatrix=DesignMatrix(:,end-length(AFNILabels)+1:end);
        useInd=ismember(AFNILabels,ConditionNames);
        DesignMatrix=array2table(UseDesignMatrix,'VariableNames',AFNILabels);
        if Compute_CondReg_CorrMat==1        
            CondReg_CorrMat=corr(UseDesignMatrix(:,useInd));
        else
            CondReg_CorrMat=[];
        end
            BpullInd=(find(useInd)*2)-1;
            TpullInd=find(useInd)*2;
        [~,BrikMat,~,~]=BrikLoad([AfniWorkDir,'glmMap+tlrc.BRIK']);
        tBrikMat=BrikMat(:,:,:,TpullInd);
        tBrikMat=reshape(tBrikMat,[size(tBrikMat,1)*size(tBrikMat,2)*size(tBrikMat,3),size(tBrikMat,4)]);
        Condition_tVals=single(tBrikMat(brain_mask(:)~=0,:));
        bBrikMat=BrikMat(:,:,:,BpullInd);
        bBrikMat=reshape(bBrikMat,[size(bBrikMat,1)*size(bBrikMat,2)*size(bBrikMat,3),size(bBrikMat,4)]);
        Condition_bVals=single(bBrikMat(brain_mask(:)~=0,:));        
    else
        numRuns=size(boldTCs,1);
        jobs = 1;
        if numRuns>1
            runStartVol=zeros(numRuns,1);
            runStartVal=0;
            for i = 1:numRuns
                  runStartVol(i,1)=runStartVal;
                  runStartVal=runStartVal+size(boldTCs{i,1},1);
            end
        else
            runStartVol=[];
            runStartVal=size(boldTCs{1,1},1);
        end 
        numVols=runStartVal;
        [ runCode ] = AFNI_3dDeconvolve(inputNames,Event_Info,[AltAfniWorkDir,'glmMap'],...
            'polort','A',...
            'TR',TR,...
            'AM2_stimInfo',TC_Info,...
            'BaseLine_stimInfo',BaseLine_Info,...
            'jobs',jobs,...
            'overwrite',1,...
            'HDR','GAM',...
            'noData',0,...
            'numVol',numVols,...
            'StartVol',runStartVol,...
            'Use3DReml',Use3DReml);
        AfniInfo.runCode=runCode;
        afniCodeFileName=[AfniWorkDir,'glmCode.tcsh'];
        ALTafniCodeFileName=[AltAfniWorkDir,'glmCode.tcsh'];
        fileID = fopen(afniCodeFileName,'w');
        fprintf(fileID,'%s',runCode);
        fclose(fileID);
        if ispc
            [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);
        else
            [~,AfniInfo.afniOut]=system(['tcsh ',ALTafniCodeFileName]);
        end

    % Load AFNI brik/head output into matlab and extract condition tVals/betas
        [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);
        numAfniRegs=size(DesignMatrix,2)-(length(AFNILabels)+size(BaseLine_Info,1));
        AfniNames=cell(numAfniRegs,1);
        for i = 1:numAfniRegs
            AfniNames{i,1}=['AFNI',num2str(i)];
        end
        if isempty(BaseLine_Info)
            AllRegressorNames=[AfniNames;AFNILabels(:)];
        else
            AllRegressorNames=[AfniNames;BaseLine_Info(:,1);AFNILabels(:)];
        end
        UseDesignMatrix=DesignMatrix(:,end-length(AFNILabels)+1:end);
        useInd=ismember(AFNILabels,ConditionNames);
        %DesignMatrix=array2table(UseDesignMatrix,'VariableNames',AFNILabels);
        if Compute_CondReg_CorrMat==1        
            CondReg_CorrMat=corr(UseDesignMatrix(:,useInd));
        else
            CondReg_CorrMat=[];
        end
        BpullInd=(find(useInd)*2)-1;
        TpullInd=find(useInd)*2;
        DesignMatrix=array2table(DesignMatrix,'VariableNames',AllRegressorNames);
        [~,BrikMat,~,~]=BrikLoad([AfniWorkDir,'glmMap+tlrc.BRIK']);
        BrikMat(:,:,:,1)=[];
        tBrikMat=BrikMat(:,:,:,TpullInd);
        tBrikMat=reshape(tBrikMat,[size(tBrikMat,1)*size(tBrikMat,2)*size(tBrikMat,3),size(tBrikMat,4)]);
        Condition_tVals=single(tBrikMat(brain_mask(:)~=0,:));
        bBrikMat=BrikMat(:,:,:,BpullInd);
        bBrikMat=reshape(bBrikMat,[size(bBrikMat,1)*size(bBrikMat,2)*size(bBrikMat,3),size(bBrikMat,4)]);
        Condition_bVals=single(bBrikMat(brain_mask(:)~=0,:)); 
    end  
    if AFNISkipDelete==0
        delete([AfniWorkDir,'*']);
    end
elseif UseAfni==2
    numRuns=size(boldTCs,1);
    jobs = 1;
    if numRuns>1
        runStartVol=zeros(numRuns,1);
        runStartVal=0;
        for i = 1:numRuns
              runStartVol(i,1)=runStartVal;
              runStartVal=runStartVal+size(boldTCs{i,1},1);
        end
    else
        runStartVol=[];
    end
    [BoldTC,brain_mask] = GLM_boldTCPrep(boldTCs,boldMasks,...
        'ParcellationMask',ParcellationMask,...
        'ResampleSizesBoldTC',ResampleSizesBoldTC,...
        'NormWithinRun',NormWithinRun,...
        'NormAcrossRun',NormAcrossRun); 
    numVols=size(BoldTC,1);
    if ispc
        if strcmp(AfniWorkDir(1,2),':')
            DirLetter=lower(AfniWorkDir(1,1));
            AltAfniWorkDir=strrep(AfniWorkDir,AfniWorkDir(1,1:2),['/mnt/',DirLetter]);
        end
    else
        AltAfniWorkDir=AfniWorkDir;
    end   
    if iscell(ConfoundTCs) && isempty(ConfoundTCs{1,1})
        ConfoundTCs=[];
    end
    if iscell(ConfoundTCsBySubj) && isempty(ConfoundTCsBySubj{1,1})
        ConfoundTCsBySubj=[];
    end        
    % Create timecourse files
    if ~isempty(ConfoundTCs) || ~isempty(ConfoundTCsBySubj)
        for run = 1:numRuns
            if ~isempty(ConfoundTCsBySubj) 
                if ~isempty(ConfoundTCsBySubj{1,1}) 
                    if size(ConfoundTCsBySubj{run,1},1)>ResampleSizes{run,1}
                        ConfoundTCsBySubj{run,1}=ConfoundTCsBySubj{run,1}(1:ResampleSizes{run,1},:);
                    end
                end
            end
            if ~isempty(ConfoundTCs) 
                if ~isempty(ConfoundTCs{1,1}) 
                    if size(ConfoundTCs{run,1},1)>ResampleSizes{run,1}
                        ConfoundTCs{run,1}=ConfoundTCs{run,1}(1:ResampleSizes{run,1},:);
                    end
                end
            end   
        end
        [tempRegressors,tempRegressorNames,~,~,tempRegressorNamesRaw] = GLM_RegressorPrep(Events,...
            TimeCourses,...
            ConfoundTCs,...
            'AM2',AM2,...
            'AM2_Events',AM2_Events,...
            'TimecourseShift',TimecourseShift,...
            'Normalize',Normalize,...
            'ConfoundTCsBySubj',ConfoundTCsBySubj,...
            'ResampleSizes',ResampleSizes);          
        AllConfoundTCNames=[];
        tempCell=cell(numRuns,1);
        for run=1:numRuns
            if ~isempty(ConfoundTCsBySubj)
                if ~isempty(ConfoundTCsBySubj{1,1})
                    AllConfoundTCNames=unique([AllConfoundTCNames;ConfoundTCsBySubj{run,2}]);
                    tempCell{run,1}=ConfoundTCsBySubj{run,1}(:,1);
                end
            end
            if ~isempty(ConfoundTCs)
                if ~isempty(ConfoundTCs{1,1})
                    AllConfoundTCNames=unique([AllConfoundTCNames;ConfoundTCs{run,2}]);
                    tempCell{run,1}=ConfoundTCs{run,1}(:,1);
                end
            end          
        end
        [~,CatInd] = Cell2CatMat(tempCell);
        ConfoundTCs=[CatMat2Cell(tempRegressors(:,ismember(tempRegressorNamesRaw,AllConfoundTCNames)),CatInd),repmat({tempRegressorNames(ismember(tempRegressorNamesRaw,AllConfoundTCNames),:)},[numRuns,1])];
        numConfoundTCs=size(ConfoundTCs{1,1},2);
        BaseLine_Info=cell(numConfoundTCs,2);
        for ConfoundTCNum = 1:numConfoundTCs
            tempConfoundTC = cell(numRuns,1);
            for run=1:numRuns
                try
                    tempConfoundTC{run,1} = ConfoundTCs{run,1}(:,ConfoundTCNum);
                catch
                    tempConfoundTC{run,1} = zeros(size(ConfoundTCs{run,1},1),1);
                end
            end
            BaseLine_Info{ConfoundTCNum,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
%             AFNILabels{AFNILabelCount,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
%             AFNILabelCount=AFNILabelCount+1;
            [BaseLine_Info{ConfoundTCNum,2}] = AFNIonsetMaker_stimfile(tempConfoundTC,AfniWorkDir,['ConfoundTC',num2str(ConfoundTCNum)],AltAfniWorkDir);
        end
    else
        BaseLine_Info=[];
    end     
    AFNILabels=cell(1);
    AFNILabelCount=1;  
    if iscell(Events)
        if isempty(Events{1,1})
            Events=[];
        end
    end
    if ~isempty(Events)
        numEvents=size(Events{1,1},2); 
        Event_Info=cell(numEvents,2);
        for eventNum = 1:numEvents
            tempEvents=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempEvents{run,1}=Events{run,1}(tInd,eventNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
            end
            Event_Info{eventNum,1}=Events{1,2}{eventNum,1};
            AFNILabels{AFNILabelCount,1}=Events{1,2}{eventNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [Event_Info{eventNum,2}] = AFNIonsetMaker_event(tempEvents,tempOnsets,AfniWorkDir,['Event',num2str(eventNum)],AltAfniWorkDir);
        end
    else
        Event_Info=[];
    end
    if iscell(TimeCourses)
        if isempty(TimeCourses{1,1})
            TimeCourses=[];
        end
    end    
    if ~isempty(TimeCourses)
        numTCs=size(TimeCourses{1,1},2);
        numAM2Events=size(AM2_Events{1,1},2);
        TC_Info=cell(numTCs,2);
        for tcNum = 1:numTCs
            tempTimeCourses=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            tempAM2_Events=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,tcNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
                if numAM2Events==0
                    tempAM2_Events{run,1}=ones(length(tempTimeCourses{run,1}),1);
                elseif numAM2Events==1
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
                else
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,tcNum);
                end
            end
            AFNILabels{AFNILabelCount,1}=['AM2_Event_',TimeCourses{1,2}{tcNum,1}];
            AFNILabelCount=AFNILabelCount+1;
            
            TC_Info{tcNum,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [TC_Info{tcNum,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['TCs',num2str(tcNum)],AltAfniWorkDir);
        end
    else
        TC_Info=[];
    end    
%     if ~isempty(TimeCourses{1,1})
%         AM2TC_Info=cell(1,2);
%         tempTimeCourses=cell(numRuns,1);
%         tempOnsets=cell(numRuns,1);
%         tempAM2_Events=cell(numRuns,1);
%         for run=1:numRuns
%             [~,tInd]=unique(TrialNum{run,1});
%             tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,:);
%             tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
%             tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
%         end
%         AM2TC_Info{1,1}=TimeCourses{1,2}{1,1};
%         AFNILabels{AFNILabelCount,1}=AM2_Events{1,2}{1,1};
%         AFNILabelCount=AFNILabelCount+1;
%         AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{1,1};
%         [AM2TC_Info{1,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['AM2_TCs'],AltAfniWorkDir);
%     else
%         AM2TC_Info=[];
%     end
    [ runCode ] = AFNI_3dDeconvolve([],Event_Info,[AltAfniWorkDir,'glmMap'],...
        'polort','A',...
        'TR',TR,...
        'AM2_stimInfo',TC_Info,...
        'BaseLine_stimInfo',BaseLine_Info,...
        'jobs',jobs,...
        'overwrite',1,...
        'HDR','GAM',...
        'noData',1,...
        'numVol',numVols,...
        'StartVol',runStartVol);
    AfniInfo.runCode=runCode;
    afniCodeFileName=[AfniWorkDir,'glmCode.tcsh'];
    ALTafniCodeFileName=[AltAfniWorkDir,'glmCode.tcsh'];
    fileID = fopen(afniCodeFileName,'w');
    fprintf(fileID,'%s',runCode);
    fclose(fileID);
    if ispc
        [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);
    else
        [~,AfniInfo.afniOut]=system(['tcsh ',ALTafniCodeFileName]);
    end
% Load AFNI brik/head output into matlab and extract condition tVals/betas
    [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);
    Regressors=DesignMatrix;
    numAfniRegs=size(DesignMatrix,2)-(length(AFNILabels)+size(BaseLine_Info,1));
    AfniNames=cell(numAfniRegs,1);
    for i = 1:numAfniRegs
        AfniNames{i,1}=['AFNI',num2str(i)];
    end
    if isempty(BaseLine_Info)
        AllRegressorNames=[AfniNames;AFNILabels(:)];
    else
        AllRegressorNames=[AfniNames;BaseLine_Info(:,1);AFNILabels(:)];
    end    
    [~,ModeFilter]=mode(DesignMatrix,1);
    EmptyRegressors=ModeFilter==size(DesignMatrix,1);
    EmptyRegressors(1,1)=0;
    RegressorNames=AllRegressorNames(EmptyRegressors==0);
    removeInd=find(sum(corr(Regressors)==1,1)>1);
    removeInd=removeInd(1,2:end);
    Regressors(:,removeInd)=[];
    AllRegressorNames(removeInd,:)=[];
    RegressorNames(removeInd,:)=[];
    EmptyRegressors(:,removeInd)=[];
    UseInd=ismember([AllRegressorNames(:);ContrastNames(:)],ConditionNames(:));     
    if strcmpi(Normalize,'zscore')
        Regressors=zscore(Regressors);
        [~,ModeFilter]=mode(Regressors,1);
        Constant=find(ModeFilter==size(Regressors,1));
        if ~isempty(Constant)
            Regressors(:,Constant(1,1))=ones(size(Regressors,1),1);
        end
    end
    if UseAfniConfounds==0
        addRegressorNames=cell(numRuns,1);
        addRegressors=zeros(size(BoldTC,1),numRuns);
        AFNIInd=contains(RegressorNames,'AFNI');
        Regressors(:,AFNIInd)=[];
        UseInd(AFNIInd,:)=[];
        UseInd=[zeros(numRuns,1);UseInd];
        EmptyRegressors(:,AFNIInd)=[];
        EmptyRegressors=[zeros(1,numRuns),EmptyRegressors];        
        RegressorNames(AFNIInd,:)=[];
        for r=1:numRuns
            addRegressorNames{r,1}=['Run',num2str(r)];
            if numRuns==1
                addRegressors=addRegressors+1;
            elseif r == numRuns
                addRegressors(runStartVol(r,1)+1:end,r)=1;
            else
                addRegressors(runStartVol(r,1)+1:runStartVol(r+1,1),r)=1;
            end
        end
        Regressors=[addRegressors,Regressors];
        RegressorNames=[addRegressorNames;RegressorNames];
    end
    if parGLM==1
        [Condition_tVals,DesignMatrix,Condition_bVals] = FastGLM_parfor(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);    
    else    
        [Condition_tVals,DesignMatrix,Condition_bVals] = FastGLM(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);
    end
    if ~isreal(Condition_tVals)
        Condition_tVals=real(Condition_tVals);
    end
    if ~isreal(Condition_bVals)
        Condition_bVals=real(Condition_bVals);
    end    
    CondReg_CorrMat=[];
    if Compute_CondReg_CorrMat==1
        CondReg_CorrMat=corrcoef(Regressors(:,UseInd==1));
    end
    DesignMatrix=array2table(DesignMatrix,'VariableNames',RegressorNames);

    Condition_tVals=Condition_tVals(:,UseInd==1); 
    Condition_bVals=Condition_bVals(:,UseInd==1); 
    delete([AfniWorkDir,'*']);
elseif UseAfni ==  3
    [MatlabRegressors,MatlabRegressorNames,~,MatlabEmptyRegressors] = GLM_RegressorPrep(Events,...
        TimeCourses,...
        ConfoundTCs,...
        'AM2',AM2,...
        'AM2_Events',AM2_Events,...
        'TimecourseShift',TimecourseShift,...
        'Normalize',Normalize,...
        'ConfoundTCsBySubj',ConfoundTCsBySubj,...
        'ResampleSizes',ResampleSizes);  
    MatlabAllRegressorNames=[MatlabRegressorNames;ContrastNames];
    MatlabUseInd=ismember(MatlabAllRegressorNames,ConditionNames);
    MatlabRegressorNames=MatlabRegressorNames(MatlabEmptyRegressors==0);    
    numRuns=size(boldTCs,1);
    jobs = 1;   
    if numRuns>1
        runStartVol=zeros(numRuns,1);
        runStartVal=0;
        for i = 1:numRuns
              runStartVol(i,1)=runStartVal;
              runStartVal=runStartVal+size(boldTCs{i,1},1);
        end
    else
        runStartVol=[];
    end
    [BoldTC,brain_mask] = GLM_boldTCPrep(boldTCs,boldMasks,...
        'ParcellationMask',ParcellationMask,...
        'ResampleSizesBoldTC',ResampleSizesBoldTC,...
        'NormWithinRun',NormWithinRun,...
        'NormAcrossRun',NormAcrossRun); 
    numVols=size(BoldTC,1);
    if ispc
        if strcmp(AfniWorkDir(1,2),':')
            DirLetter=lower(AfniWorkDir(1,1));
            AltAfniWorkDir=strrep(AfniWorkDir,AfniWorkDir(1,1:2),['/mnt/',DirLetter]);
        end
    else
        AltAfniWorkDir=AfniWorkDir;
    end
    AFNILabels=cell(1);
    AFNILabelCount=1;    
    if ~isempty(Events)
        numEvents=size(Events{1,1},2); 
        Event_Info=cell(numEvents,2);
        for eventNum = 1:numEvents
            tempEvents=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempEvents{run,1}=Events{run,1}(tInd,eventNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
            end
            Event_Info{eventNum,1}=Events{1,2}{eventNum,1};
            AFNILabels{AFNILabelCount,1}=Events{1,2}{eventNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [Event_Info{eventNum,2}] = AFNIonsetMaker_event(tempEvents,tempOnsets,AfniWorkDir,['Event',num2str(eventNum)],AltAfniWorkDir);
        end
    else
        Event_Info=[];
    end
    if ~isempty(TimeCourses{1,1})
        numTCs=size(TimeCourses{1,1},2);
        numAM2Events=size(AM2_Events{1,1},2);
        TC_Info=cell(numTCs,2);
        for tcNum = 1:numTCs
            tempTimeCourses=cell(numRuns,1);
            tempOnsets=cell(numRuns,1);
            tempAM2_Events=cell(numRuns,1);
            for run=1:numRuns
                [~,tInd]=unique(TrialNum{run,1});
                tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,tcNum);
                tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
                if numAM2Events==0
                    tempAM2_Events{run,1}=ones(length(tempTimeCourses{run,1}),1);
                elseif numAM2Events==1
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
                else
                    tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,tcNum);
                end
            end
            AFNILabels{AFNILabelCount,1}=['AM2_Event_',TimeCourses{1,2}{tcNum,1}];
            AFNILabelCount=AFNILabelCount+1;
            
            TC_Info{tcNum,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{tcNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [TC_Info{tcNum,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['TCs',num2str(tcNum)],AltAfniWorkDir);
        end
    else
        TC_Info=[];
    end    
    [ runCode ] = AFNI_3dDeconvolve([],Event_Info,[AltAfniWorkDir,'glmMap'],...
        'polort','A',...
        'TR',TR,...
        'AM2_stimInfo',TC_Info,...
        'BaseLine_stimInfo',BaseLine_Info,...
        'jobs',jobs,...
        'overwrite',1,...
        'HDR','GAM',...
        'noData',1,...
        'numVol',numVols,...
        'StartVol',runStartVol);
    AfniInfo.runCode=runCode;
    afniCodeFileName=[AfniWorkDir,'glmCode.tcsh'];
    ALTafniCodeFileName=[AltAfniWorkDir,'glmCode.tcsh'];
    fileID = fopen(afniCodeFileName,'w');
    fprintf(fileID,'%s',runCode);
    fclose(fileID);
    if ispc
        [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);
    else
        [~,AfniInfo.afniOut]=system(['tcsh ',ALTafniCodeFileName]);
    end
% Load AFNI brik/head output into matlab and extract condition tVals/betas
    [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);
    Regressors=DesignMatrix;
    numAfniRegs=size(DesignMatrix,2)-(length(AFNILabels)+size(BaseLine_Info,1));
    AfniNames=cell(numAfniRegs,1);
    for i = 1:numAfniRegs
        AfniNames{i,1}=['AFNI',num2str(i)];
    end
    AllRegressorNames=[AfniNames;BaseLine_Info(:,1);AFNILabels(:)];
    [~,ModeFilter]=mode(DesignMatrix,1);
    EmptyRegressors=ModeFilter==size(DesignMatrix,1);
    EmptyRegressors(1,1)=0;
    RegressorNames=AllRegressorNames(EmptyRegressors==0);
    removeInd=find(sum(corr(Regressors)==1,1)>1);
    removeInd=removeInd(1,2:end);
    Regressors(:,removeInd)=[];
    AllRegressorNames(removeInd,:)=[];
    RegressorNames(removeInd,:)=[];
    EmptyRegressors(:,removeInd)=[];
    UseInd=ismember([AllRegressorNames(:);ContrastNames(:)],ConditionNames(:));     
    if strcmpi(Normalize,'zscore')
        Regressors=zscore(Regressors);
        [~,ModeFilter]=mode(Regressors,1);
        Constant=find(ModeFilter==size(Regressors,1));
        if ~isempty(Constant)
            Regressors(:,Constant(1,1))=ones(size(Regressors,1),1);
        end
    end

    if parGLM==1
        [Condition_tVals,DesignMatrix] = FastGLM_parfor(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);    
    else    
        [Condition_tVals,DesignMatrix] = FastGLM(BoldTC,Regressors,Contrasts,BatchSize,...
            'EmptyRegressors',EmptyRegressors,...
            'ResampleToFit',ResampleToFit,...
            'UseStat',UseStat);
    end
    if ~isreal(Condition_tVals)
        Condition_tVals=real(Condition_tVals);
    end
    CondReg_CorrMat=[];
    if Compute_CondReg_CorrMat==1
        CondReg_CorrMat=corrcoef(Regressors(:,UseInd==1));
    end
    DesignMatrix=array2table(DesignMatrix,'VariableNames',RegressorNames);

    Condition_tVals=Condition_tVals(:,UseInd==1);     
    delete([AfniWorkDir,'*']);   
    
end
end