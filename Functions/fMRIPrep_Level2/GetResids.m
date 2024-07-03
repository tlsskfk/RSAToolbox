function [Resids,brain_mask,AllRegressorNames,DesignMatrix,AfniInfo] = GetResids(Events,TimeCourses,ConfoundTCs,boldTCs,boldMasks,varargin)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AfniInfo=[];
if UseAfni==0
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

    RegressorNames=RegressorNames(EmptyRegressors==0);
    [Resids] = FastOLSRegress_Resids(BoldTC,Regressors);    
    DesignMatrix=array2table(Regressors,'VariableNames',AllRegressorNames);    
elseif UseAfni==1
    numRuns=size(boldTCs,1);
    jobs = feature('numcores');
    if strcmp(AfniWorkDir(1,2),':')
        DirLetter=lower(AfniWorkDir(1,1));
        AltAfniWorkDir=strrep(AfniWorkDir,AfniWorkDir(1,1:2),['/mnt/',DirLetter]);
    end    
    AFNILabels=cell(1);
    AFNILabelCount=1; 
    % Create timecourse files
    if ~isempty(ConfoundTCs)
        numConfoundTCs=size(ConfoundTCs{1,1},2);  
        BaseLine_Info=cell(numConfoundTCs,2);
        for ConfoundTCNum = 1:numConfoundTCs
            tempConfoundTC = cell(numRuns,1);
            for run=1:numRuns
                tempConfoundTC{run,1} = ConfoundTCs{run,1}(:,ConfoundTCNum);
            end
            BaseLine_Info{ConfoundTCNum,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
            AFNILabels{AFNILabelCount,1} = ConfoundTCs{1,2}{ConfoundTCNum,1};
            AFNILabelCount=AFNILabelCount+1;
            [BaseLine_Info{ConfoundTCNum,2}] = AFNIonsetMaker_stimfile(tempConfoundTC,AfniWorkDir,['ConfoundTC',num2str(ConfoundTCNum)],AltAfniWorkDir);
        end
    else
        BaseLine_Info=[];
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
    if ~isempty(TimeCourses{1,1})
        AM2TC_Info=cell(1,2);
        tempTimeCourses=cell(numRuns,1);
        tempOnsets=cell(numRuns,1);
        tempAM2_Events=cell(numRuns,1);
        for run=1:numRuns
            [~,tInd]=unique(TrialNum{run,1});
            tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,:);
            tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
            tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
        end
        AM2TC_Info{1,1}=TimeCourses{1,2}{1,1};
        AFNILabels{AFNILabelCount,1}=AM2_Events{1,2}{1,1};
        AFNILabelCount=AFNILabelCount+1;
        AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{1,1};
        [AM2TC_Info{1,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['AM2_TCs'],AltAfniWorkDir);
    else
        AM2TC_Info=[];
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
            'GOFORIT',12,...
            'AM2_stimInfo',AM2TC_Info,...
            'BaseLine_stimInfo',BaseLine_Info,...
            'residOut',[AltAfniWorkDir,'Resids'],...
            'jobs',jobs,...
            'overwrite',1,...
            'HDR','GAM');
        AfniInfo.runCode=runCode;
        afniCodeFileName=[AfniWorkDir,'glmCode.tcsh'];
        ALTafniCodeFileName=[AltAfniWorkDir,'glmCode.tcsh'];
        fileID = fopen(afniCodeFileName,'w');
        fprintf(fileID,'%s',runCode);
        fclose(fileID);
        [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);

    % Load AFNI brik/head output into matlab and extract condition tVals/betas
        [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);

        UseDesignMatrix=DesignMatrix(:,end-length(AFNILabels)+1:end);
        DesignMatrix=array2table(UseDesignMatrix,'VariableNames',AFNILabels);

        [~,BrikMat,~,~]=BrikLoad([AfniWorkDir,'Resids+tlrc.BRIK']);
        BrikMat=reshape(BrikMat,[size(BrikMat,1)*size(BrikMat,2)*size(BrikMat,3),size(BrikMat,4)]);
        BrikMat=BrikMat';
        tempMask=brain_mask(:)';
        Resids=BrikMat(:,find(tempMask));
        brain_mask(brain_mask==1)=[1:sum(brain_mask(:))]; 
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
            'AM2_stimInfo',AM2TC_Info,...
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
        [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);

    % Load AFNI brik/head output into matlab and extract condition tVals/betas
        [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);
        numAfniRegs=size(DesignMatrix,2)-(length(AFNILabels)+size(BaseLine_Info,1));
        AfniNames=cell(numAfniRegs,1);
        for i = 1:numAfniRegs
            AfniNames{i,1}=['AFNI',num2str(i)];
        end
        AllRegressorNames=[AfniNames;BaseLine_Info(:,1);AFNILabels(:)];
        UseDesignMatrix=DesignMatrix(:,end-length(AFNILabels)+1:end);
        useInd=ismember(AFNILabels,ConditionNames);
        %DesignMatrix=array2table(UseDesignMatrix,'VariableNames',AFNILabels);
        if Compute_CondReg_CorrMat==1        
            CondReg_CorrMat=corr(UseDesignMatrix(:,useInd));
        else
            CondReg_CorrMat=[];
        end
        if strcmpi(UseStat,'B')
            pullInd=(find(useInd)*2)-1;
        else
            pullInd=find(useInd)*2;
        end
        DesignMatrix=array2table(DesignMatrix,'VariableNames',AllRegressorNames);
        [~,BrikMat,~,~]=BrikLoad([AfniWorkDir,'glmMap+tlrc.BRIK']);
        BrikMat(:,:,:,1)=[];
        BrikMat=BrikMat(:,:,:,pullInd);
        BrikMat=reshape(BrikMat,[size(BrikMat,1)*size(BrikMat,2)*size(BrikMat,3),size(BrikMat,4)]);
        Condition_tVals=single(BrikMat(brain_mask(:)~=0,:));
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
    if strcmp(AfniWorkDir(1,2),':')    
        DirLetter=lower(AfniWorkDir(1,1));
        AltAfniWorkDir=strrep(AfniWorkDir,AfniWorkDir(1,1:2),['/mnt/',DirLetter]);
    end    

    % Create timecourse files
    if ~isempty(ConfoundTCs)
        numConfoundTCs=size(ConfoundTCs{1,1},2);  
        BaseLine_Info=cell(numConfoundTCs,2);
        for ConfoundTCNum = 1:numConfoundTCs
            tempConfoundTC = cell(numRuns,1);
            for run=1:numRuns
                tempConfoundTC{run,1} = ConfoundTCs{run,1}(:,ConfoundTCNum);
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
    if ~isempty(TimeCourses{1,1})
        AM2TC_Info=cell(1,2);
        tempTimeCourses=cell(numRuns,1);
        tempOnsets=cell(numRuns,1);
        tempAM2_Events=cell(numRuns,1);
        for run=1:numRuns
            [~,tInd]=unique(TrialNum{run,1});
            tempTimeCourses{run,1}=TimeCourses{run,1}(tInd,:);
            tempOnsets{run,1}=TrialOnsetTimes{run,1}(tInd,1);
            tempAM2_Events{run,1}=AM2_Events{run,1}(tInd,1);
        end
        AM2TC_Info{1,1}=TimeCourses{1,2}{1,1};
        AFNILabels{AFNILabelCount,1}=AM2_Events{1,2}{1,1};
        AFNILabelCount=AFNILabelCount+1;
        AFNILabels{AFNILabelCount,1}=TimeCourses{1,2}{1,1};
        [AM2TC_Info{1,2}] = AFNIonsetMaker_AM2TC(tempTimeCourses,tempAM2_Events,tempOnsets,AfniWorkDir,['AM2_TCs'],AltAfniWorkDir);
    else
        AM2TC_Info=[];
    end
    [ runCode ] = AFNI_3dDeconvolve([],Event_Info,[AltAfniWorkDir,'glmMap'],...
        'polort','A',...
        'TR',TR,...
        'AM2_stimInfo',AM2TC_Info,...
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
    [~,AfniInfo.afniOut]=system(['wsl tcsh ',ALTafniCodeFileName]);

% Load AFNI brik/head output into matlab and extract condition tVals/betas
    [~,DesignMatrix, ~] = Read_1D ([AfniWorkDir,'glmMap.xmat']);
    Regressors=DesignMatrix;
    if isempty(AFNILabels{1,1})
        numAfniRegs=size(DesignMatrix,2)-size(BaseLine_Info,1);
    else
        numAfniRegs=size(DesignMatrix,2)-(length(AFNILabels)+size(BaseLine_Info,1));
    end
    AfniNames=cell(numAfniRegs,1);
    for i = 1:numAfniRegs
        AfniNames{i,1}=['AFNI',num2str(i)];
    end
    if isempty(AFNILabels{1,1})
        AllRegressorNames=[AfniNames;BaseLine_Info(:,1)];
    else
        AllRegressorNames=[AfniNames;BaseLine_Info(:,1);AFNILabels(:)];
    end
    [~,ModeFilter]=mode(DesignMatrix,1);
    EmptyRegressors=ModeFilter==size(DesignMatrix,1);
    EmptyRegressors(1,1)=0;
 
    if strcmpi(Normalize,'zscore')
        Regressors=zscore(Regressors);
        [~,ModeFilter]=mode(Regressors,1);
        Constant=find(ModeFilter==size(Regressors,1));
        if ~isempty(Constant)
            Regressors(:,Constant(1,1))=ones(size(Regressors,1),1);
        end
    end

    [Resids] = FastOLSRegress_Resids(BoldTC,Regressors);    
    DesignMatrix=array2table(Regressors,'VariableNames',AllRegressorNames);    
    
    delete([AfniWorkDir,'*']);
end

end
