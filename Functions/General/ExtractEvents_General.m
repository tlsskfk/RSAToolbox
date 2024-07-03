function [ TrialData,fMRI_Events,InputParameters ] = ExtractEvents_General( Stimuli,Onsets,varargin)
%Written by David Rothlein
% Convert behavioral or timecourse data into multiple usable formats.
% Input: 
    %Stimuli - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %Onsets - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.
        
% Output: 
    %TrialData - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %fMRI_Events - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.        
    %InputParameters - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus. 
        
% Optional input
    %condNames - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %ExemplarTypes - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.
    %Responses - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %respTypes - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.
    %EventTypes - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %BlockThresh - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.
    %AddTimeConstant - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %mseq - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.   
    %fMRIDur - a T by 1 cell or matrix where T is the number of trials. 
        %Contains strings or values indicating different stimuli.
    %TimeCourses - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.  
    %TimeIncrementInSec - a T by 1 cell or matrix where T is the number of trials.
        %Contains the onset time (in seconds) of each stimulus.        
        
%If input Stimuli is a vector, convert to cell of same size        
if ~iscell(Stimuli)
    Stimuli=num2cell(Stimuli);
end

%Make sure Stimuli and onsets are vertical
Stimuli=Stimuli(:);
Onsets=Onsets(:);

FilterInd=isnan(Onsets); %Identify and remove nans from input
%if optional input CondList is not set, make it the same as Stimuli
CondList = VariableSetter( 'CondList',Stimuli,varargin); 
%Remove NANs
Stimuli(FilterInd,:)=[]; 
CondList(FilterInd,:)=[];

%If optional variables are not set, set as defaults.
TrialDurations = VariableSetter( 'TrialDurations',[],varargin);
condNames = VariableSetter( 'condNames',[],varargin);
ExemplarTypes = VariableSetter( 'ExemplarTypes',[],varargin);
Responses = VariableSetter( 'Responses',[],varargin);
respTypes = VariableSetter( 'respTypes',[],varargin);
EventTypes = VariableSetter( 'EventTypes',[],varargin);
BlockThresh = VariableSetter( 'BlockThresh',[],varargin);
AddTimeConstant = VariableSetter( 'AddTimeConstant',[],varargin);
mseq = VariableSetter( 'mseq',[0],varargin);
fMRIDur = VariableSetter( 'fMRIDur',[],varargin);
TimeCourses = VariableSetter( 'TimeCourses',[],varargin);
TimeIncrementInSec = VariableSetter( 'TimeIncrementInSec',0.5,varargin);
UseParameters = VariableSetter( 'UseParameters',[],varargin);
SetUseParameters=VariableSetter( 'SetUseParameters',0,varargin);
%% Define optional inputs using GUI
if isempty(condNames) && isempty(UseParameters)
    [~,condLabels,condGroups] = uiBehaviorTrialLabels(CondList,'Assign condition labels to trial types');
    condNames={condLabels,condGroups};
elseif isempty(UseParameters)
    condLabels=condNames{1,1};
    condGroups=condNames{1,2};
else
    condNames=UseParameters.condNames;
    condLabels=condNames{1,1};
    condGroups=condNames{1,2};   
end
numConds=length(condLabels);

if isempty(ExemplarTypes) && isempty(UseParameters)
    [~,ExemplarLabels,ExemplarGroups] = uiBehaviorTrialLabels(Stimuli,'Assign exemplar labels to trial types');
    ExemplarTypes={ExemplarLabels,ExemplarGroups};
elseif isempty(UseParameters)
    ExemplarLabels=ExemplarTypes{1,1};
    ExemplarGroups=ExemplarTypes{1,2};
else
    ExemplarTypes=UseParameters.ExemplarTypes;
    ExemplarLabels=ExemplarTypes{1,1};
    ExemplarGroups=ExemplarTypes{1,2};
end
numExemplars=length(ExemplarTypes{1,1});

if ~isempty(Responses)
    if isempty(respTypes) && isempty(UseParameters)
        [~,respLabels,respGroups] = uiBehaviorTrialLabels(Responses,'Assign labels to recorded responses');
        respTypes={respLabels,respGroups};
    elseif isempty(UseParameters)
        respLabels=respTypes{1,1};
        respGroups=respTypes{1,2};
    else
        respTypes=UseParameters.respTypes;
        respLabels=respTypes{1,1};
        respGroups=respTypes{1,2};        
    end   
    if isempty(EventTypes) && isempty(UseParameters)
        [~,EventLabels,EventGroups] = uiBehaviorTrialLabels([condLabels,respLabels],'Assign event labels to condition/response pairs');
        EventTypes={EventLabels,EventGroups};
    elseif isempty(UseParameters)
        EventLabels=EventTypes{1,1};
        EventGroups=EventTypes{1,2};
    else
        EventTypes=UseParameters.EventTypes;
        EventLabels=EventTypes{1,1};
        EventGroups=EventTypes{1,2};
    end
    numEvents=length(EventLabels);
    numResp=length(respLabels);    
end
if isempty(UseParameters)
    if isempty(BlockThresh)
        [BlockThresh] = uiEnterName('8','Enter trial duration for block designation in seconds');
        BlockThresh = str2num(BlockThresh);
    end
    if isempty(AddTimeConstant)
        [AddTimeConstant] = uiEnterName('0','Enter time duration to add to all onsets');
        AddTimeConstant = str2num(AddTimeConstant);
    end
else
    BlockThresh=UseParameters.BlockThresh;
    AddTimeConstant=UseParameters.AddTimeConstant;
end
if SetUseParameters==0
    %% Identify trial and experimental durations
    numTrials=length(Stimuli);
    TrialData.TrialNum=[1:numTrials]';
    TrialData.TrialOnset=Onsets+AddTimeConstant; 

    %Duration will not be defined for last trial using diff fcn. Therefore the
    %duration of the last trial will be the mean of all the other durations.
    %%%%% THIS CAN CHANGE IF A BETTER IDEA COMES ALONG %%%%%
    if isempty(TrialDurations)
        TrialData.Duration=diff(Onsets); 
        TrialData.Duration=[TrialData.Duration;mean(TrialData.Duration)];
    else
        TrialData.Duration=TrialDurations;
    end
    TrialData.ExperimentDuration=Onsets(end)+TrialData.Duration(end);
    TrialData.StartTimePad=TrialData.Duration;
    if ~isempty(fMRIDur)
        TrialData.EndTimePad=fMRIDur-TrialData.ExperimentDuration; 
    else
        TrialData.EndTimePad=0;
    end

    %% Define basic trial types and create the Trial2Time Vector
    %Create matrix to convert trial space into time space at a 100ms
    %resolution
    TrialData.Trial2Time=zeros(round(TrialData.ExperimentDuration*10),1);

    TrialData.Block=single(TrialData.Duration>=BlockThresh);
    TrialData.Event=single(TrialData.Duration<BlockThresh);
    TrialData.AverageEventDur=nanmean(TrialData.Duration(TrialData.Event==1,:),1);
    TrialData.Stimuli=Stimuli;

    Trial2TimeIndex=round(Onsets*10)+1;
    Trial2TimeIndex=[Trial2TimeIndex,[Trial2TimeIndex(2:end)-1;Trial2TimeIndex(end,1)+round(TrialData.Duration(end,1)*10)-1]];

    for trial=1:numTrials
        TrialData.Trial2Time(Trial2TimeIndex(trial,1):Trial2TimeIndex(trial,2),1)=trial;
    end

    %% Construct Condition, Exemplar, and if exists, Response and Event Vectors
    TrialData.EventLists.Condition=repmat({'null'},[numTrials,1]);
    TrialData.EventLists.Exemplar=repmat({'null'},[numTrials,1]);
    TrialData.EventLists.Response=repmat({'null'},[numTrials,1]);
    TrialData.EventLists.Events=repmat({'null'},[numTrials,1]);
    %% Create condition lists
    for condi=1:numConds
        condName=condLabels{1,condi};
        condGroup=condGroups{1,condi};
        condInd=ismember(Stimuli,condGroup);
        TrialData.EventSequences.(condName)=single(ismember(Stimuli,condGroup));
        TrialData.EventLists.Condition(condInd,1)={condName};
        TrialData.EventOnsets.(condName)=TrialData.TrialOnset(condInd,1);
        TrialData.EventDurations.(condName)=TrialData.Duration(condInd,1);
        TrialData.fMRI_Events.(condName)=trial2time(TrialData.EventSequences.(condName),TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
    end
    %%Create exemplar lists
    for condi=1:numExemplars
        condName=ExemplarLabels{1,condi};
        condGroup=ExemplarGroups{1,condi};
        condInd=ismember(Stimuli,condGroup);
        TrialData.EventSequences.(condName)=single(ismember(Stimuli,condGroup));
        TrialData.EventLists.Exemplar(condInd,1)={condName};
        TrialData.EventOnsets.(condName)=TrialData.TrialOnset(condInd,1);
        TrialData.EventDurations.(condName)=TrialData.Duration(condInd,1);
        TrialData.fMRI_Events.(condName)=trial2time(TrialData.EventSequences.(condName),TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
    end
    %% create exemplars for sequential trial pairs
    if mseq==2    
        for condi1=1:numExemplars
            condName1=ExemplarLabels{1,condi1};
            condGroup1=ExemplarGroups{1,condi1};
            condInd1=single(ismember(Stimuli,condGroup1));
            condInd1=[0;condInd1(1:end-1,1)];       
            for condi2=1:numExemplars
                condName2=ExemplarLabels{1,condi2};
                condGroup2=ExemplarGroups{1,condi2};
                condInd2=single(ismember(Stimuli,condGroup2));  
                condName=['mSeq_',condName1,'_',condName2];
                condInd=(condInd1+condInd2)==2;
                TrialData.mseq.EventSequences.(condName)=single(condInd);
                TrialData.mseq.EventLists.Exemplar(condInd,1)={condName};
                TrialData.mseq.EventOnsets.(condName)=TrialData.TrialOnset(condInd,1);
                TrialData.mseq.EventDurations.(condName)=TrialData.Duration(condInd,1);
                TrialData.mseq.fMRI_Events.(condName)=trial2time(TrialData.mseq.EventSequences.(condName),TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
            end
        end
    else
        TrialData.mseq=[];
    end
    %% create lists for responses and response-based events
    if ~isempty(Responses)
        for condi=1:numResp
            condName=respLabels{1,condi};
            condGroup=respGroups{1,condi};
            condInd=ismember(Responses,condGroup);
            TrialData.EventSequences.(condName)=single(ismember(Responses,condGroup));
            TrialData.EventLists.Response(condInd,1)={condName};
            TrialData.EventOnsets.(condName)=TrialData.TrialOnset(condInd,1);
            TrialData.EventDurations.(condName)=TrialData.Duration(condInd,1);
            TrialData.fMRI_Events.(condName)=trial2time(TrialData.EventSequences.(condName),TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
        end
        for condi=1:numEvents
            condName=EventLabels{1,condi};
            condGroup=EventGroups{1,condi};
            condInd=ismember(TrialData.EventLists.Condition,condGroup) & ismember(TrialData.EventLists.Response,condGroup);
            TrialData.EventSequences.(condName)=single(condInd);
            TrialData.EventLists.Events(condInd,1)={condName};
            TrialData.EventOnsets.(condName)=TrialData.TrialOnset(condInd,1);
            TrialData.EventDurations.(condName)=TrialData.Duration(condInd,1);
            TrialData.fMRI_Events.(condName)=trial2time(TrialData.EventSequences.(condName),TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
        end
    end
    TrialData.fMRI_Events.TrialNum=trial2time(TrialData.TrialNum,TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
    TrialData.fMRI_Events.TrialOnsetTime=trial2time(TrialData.TrialOnset,TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
    %% create output table of temporal time series conditions
    if mseq==2
        fMRI_Events=[struct2table(TrialData.fMRI_Events),struct2table(TrialData.mseq.fMRI_Events)];
    else
        fMRI_Events=struct2table(TrialData.fMRI_Events);
    end
    %% Add RTs and other timecourses
    if ~isempty(TimeCourses)
        TimeCoursesMat=TimeCourses{1,1};
        TimeCoursesLabels=TimeCourses{1,2};
        for iTC = 1:length(TimeCoursesLabels)
            TimeCourse=TimeCoursesMat(iTC,:);
            TimeCourseLabel=TimeCoursesLabels{1,iTC};
            TrialData.TimeCourses.(TimeCourseLabel)=TimeCourse;
            TrialData.fMRI_Events.(TimeCourseLabel)=trial2time(TimeCourse,TrialData.Trial2Time,fMRIDur,TimeIncrementInSec);
        end
    end
else
    TrialData=[];
    fMRI_Events=[];
end

%% compile inputparameters for additional runs
InputParameters.condNames=condNames;
InputParameters.ExemplarTypes = ExemplarTypes;
InputParameters.respTypes = respTypes;
InputParameters.EventTypes = EventTypes;
InputParameters.BlockThresh = BlockThresh;
InputParameters.AddTimeConstant = AddTimeConstant;

end

%% Convert trial sequence to time series
function [EventsMat_ReSamp]=trial2time(TrialMat,Trial2Time,ScanDurInSec,TimeIncrementInSec)
numCond=size(TrialMat,2);

%convert TrialMat from trial to decisecond space
NullTrialIndex=Trial2Time==0;
Trial2Time(NullTrialIndex)=max(Trial2Time)+1; %Remove 0 trial index and replace last trial + 1 
Trial2TimeMat=zeros(length(Trial2Time),numCond);

for cond = 1:numCond %Iterate through conditions and assign trial values in time space
    tempTrial=[TrialMat(:,cond);nan]; %pull condition trial TC & add nan to end.
                                      %nan will be assigned to former 0-trial index
    Trial2TimeMat(:,cond)=tempTrial(Trial2Time); %reindex trial values in time space
end
scanDurInMSec=round(ScanDurInSec*10);
if size(Trial2TimeMat,1)>scanDurInMSec
    Trial2TimeMat(scanDurInMSec+1:end,:)=[];
elseif size(Trial2TimeMat,1)<scanDurInMSec
    Trial2TimeMat=[Trial2TimeMat;nan(scanDurInMSec-size(Trial2TimeMat,1),size(Trial2TimeMat,2))];
end
TimeIncrementInMSec=round(TimeIncrementInSec*10);

ResampleLength=size(Trial2TimeMat,1)/TimeIncrementInMSec;
EventsMat_ReSamp = imresize(Trial2TimeMat,[ResampleLength,size(Trial2TimeMat,2)],'nearest');

end    
