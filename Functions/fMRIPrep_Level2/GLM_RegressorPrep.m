function [X,XNames,XTypes,EmptyRegressor,XNamesRaw,ConfoundTCs] = GLM_RegressorPrep(Events,TimeCourses,ConfoundTCs,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Normalize] = VariableSetter('Normalize','zscore',varargin);
[ConfoundTCsBySubj] = VariableSetter('ConfoundTCsBySubj',{[]},varargin);
[ResampleSizes] = VariableSetter('ResampleSizes',[],varargin);
[AM2]=VariableSetter('AM2',1,varargin);
[AM2_Events]=VariableSetter('AM2_Events',{[]},varargin);
TimecourseShift = VariableSetter('TimecourseShift',0,varargin);

runTC=[];
runTCNames=[];
if ~iscell(Events) && ~isempty(Events)
    Events={Events};
end
if ~iscell(TimeCourses) && ~isempty(TimeCourses)
    TimeCourses={TimeCourses};
end
if ~iscell(ConfoundTCs) && ~isempty(ConfoundTCs)
    ConfoundTCs={ConfoundTCs};
end
if ~isempty(Events)
    numRuns=size(Events,1);
elseif ~isempty(TimeCourses)
    numRuns=size(TimeCourses,1);
else
    numRuns=size(ConfoundTCs,1);
end

RunDurs=cell(numRuns,1);
if ~isempty(Events{1,1})
    for run = 1:numRuns
        Events{run,1}(isnan(Events{run,1}))=0;
        Events{run,1}=HRFConvolve(Events{run,1});
        RunDurs{run,1}=size(Events{run,1},1);
        if ~isempty(ResampleSizes)
            Events{run,1}=imresize(Events{run,1},[ResampleSizes{run,1},size(Events{run,1},2)]);
        end
    end
    if numRuns>1
        [Events,EventNames,runTC,runTCNames]=NameAlignVar(Events(:,1),Events(:,2));
    else
        EventNames=Events{1,2};
        Events=Events{1,1};
        runTC=[];
        runTCNames=[];
    end
else
    Events=[];
    EventNames=[];
end

if ~isempty(TimeCourses{1,1})
    for run = 1:numRuns
        if TimecourseShift~=0
            TimeCourses{run,1}(end-TimecourseShift+1:end,:)=[];
            TimeCourses{run,1}=[nan(TimecourseShift,size(TimeCourses{run,1},2));TimeCourses{run,1}];
            TimeCourses{run,1}=fillmissing(TimeCourses{run,1},'nearest');
        end
        if AM2~=0
            if isempty(AM2_Events{1,1})
                for TCNum=1:size(TimeCourses{run,1},2)
                    TempTC=TimeCourses{run,1}(:,TCNum);
                    AM2_Event=isnan(TempTC)==0;
                    TempTC(AM2_Event==0)=0;
                    TempTC(TempTC~=0)=TempTC(TempTC~=0)-nanmean(TempTC(TempTC~=0),1);
                    TempTC(isnan(TempTC))=0;
                    TimeCourses{run,1}(:,TCNum)=HRFConvolve(TempTC);
                end
            else
                AM2_Events{run,1}(isnan(AM2_Events{run,1}))=0;
                for TCNum=1:size(TimeCourses{run,1},2)
                    TempTC=TimeCourses{run,1}(:,TCNum);
                    try
                        AM2_Event=single(isnan(TempTC)==0+AM2_Events{run,1}(:,TCNum)==0)==0;
                    catch
                        AM2_Event=single(isnan(TempTC)==0+AM2_Events{run,1}(:,1)==0)==0;
                    end
                    TempTC(AM2_Event==0)=0;
                    TempTC(TempTC~=0)=TempTC(TempTC~=0)-nanmean(TempTC(TempTC~=0),1);
                    TempTC(isnan(TempTC))=0;
                    TimeCourses{run,1}(:,TCNum)=HRFConvolve(TempTC);
                end  
                [AM2_Events{run,2},AM2_Events_Ind]=unique(AM2_Events{run,2});
                TimeCourses{run,1}=[TimeCourses{run,1},HRFConvolve(AM2_Events{run,1}(:,AM2_Events_Ind))];
                TimeCourses{run,2}=[TimeCourses{run,2};AM2_Events{run,2}];
            end
        end
        if ~isempty(ResampleSizes)
            TimeCourses{run,1}=imresize(TimeCourses{run,1},[ResampleSizes{run,1},size(TimeCourses{run,1},2)]);
        end
        if isempty(Events)
            RunDurs{run,1}=size(TimeCourses{run,1},1);
        end
    end
    if numRuns>1
        [TimeCourses,TimeCoursesNames,runTC,runTCNames]=NameAlignVar(TimeCourses(:,1),TimeCourses(:,2));
    else
        TimeCoursesNames=TimeCourses{1,2};
        TimeCourses=TimeCourses{1,1};
        runTC=[];
        runTCNames=[];
    end
else
    TimeCourses=[];
    TimeCoursesNames=[];
end

if ~isempty(ConfoundTCs)
    if any(~cellfun(@isempty,ConfoundTCs(:,1)))

        for run = 1:numRuns
            if isempty(ConfoundTCs{run,1})
                useInd=find(~cellfun(@isempty,ConfoundTCs(:,1)));
                ConfoundTCs{run,1}=zeros(ResampleSizes{run,1},size(ConfoundTCs{useInd(1,1),1},2));
                ConfoundTCs{run,2}=ConfoundTCs{useInd(1,1),2};
            end
            ConfoundTCs{run,1}=fillmissing(ConfoundTCs{run,1},'nearest');
            if ~isempty(ResampleSizes)
                ConfoundTCs{run,1}=imresize(ConfoundTCs{run,1},[ResampleSizes{run,1},size(ConfoundTCs{run,1},2)]);
            elseif ~isempty(RunDurs{1,1}) 
                ConfoundTCs{run,1}=imresize(ConfoundTCs{run,1},[RunDurs{run,1},size(ConfoundTCs{run,1},2)],'nearest');
            end
        end
        if numRuns>1
            [~,~,runTC,runTCNames]=NameAlignVar(ConfoundTCs(:,1),ConfoundTCs(:,2));
            [ConfoundTCs,ConfoundTCsNames,ConfoundTCsNamesRaw]=cell2ConfoundMat(ConfoundTCs,Normalize);          
        else
            ConfoundTCsNames=ConfoundTCs{1,2}; 
            ConfoundTCsNamesRaw=ConfoundTCsNames;
            ConfoundTCs=ConfoundTCs{1,1};
            runTC=[];
            runTCNames=[];
        end
    else
        ConfoundTCs=[];
        ConfoundTCsNames=[];
        ConfoundTCsNamesRaw=[];
    end
else
    ConfoundTCsNames=[];
    ConfoundTCsNamesRaw=[];    
end

if iscell(ConfoundTCsBySubj)
    if ~isempty(ConfoundTCsBySubj{1,1})
        for run = 1:numRuns
            ConfoundTCsBySubj{run,1}=fillmissing(ConfoundTCsBySubj{run,1},'nearest');
            if ~isempty(ResampleSizes)
                ConfoundTCsBySubj{run,1}=imresize(ConfoundTCsBySubj{run,1},[ResampleSizes{run,1},size(ConfoundTCsBySubj{run,1},2)]);
            elseif ~isempty(RunDurs{1,1}) 
                ConfoundTCsBySubj{run,1}=imresize(ConfoundTCsBySubj{run,1},[RunDurs{run,1},size(ConfoundTCsBySubj{run,1},2)],'nearest');
            end
        end
        if numRuns>1
            [ConfoundTCsBySubj,ConfoundTCsBySubjNames,runTC,runTCNames]=NameAlignVar(ConfoundTCsBySubj(:,1),ConfoundTCsBySubj(:,2));         
        else
            ConfoundTCsBySubjNames=ConfoundTCsBySubj{1,2};        
            ConfoundTCsBySubj=ConfoundTCsBySubj{1,1};
            runTC=[];
            runTCNames=[];
        end
    else
        ConfoundTCsBySubj=[];
        ConfoundTCsBySubjNames=[];
    end
else
    ConfoundTCsBySubj=[];
    ConfoundTCsBySubjNames=[];
end
if strcmpi(Normalize,'zscore')
    ConfoundTCsBySubj=zscore(ConfoundTCsBySubj);
end
ConfoundTCs=[ConfoundTCsBySubj,ConfoundTCs];
ConfoundTCsNames=[ConfoundTCsBySubjNames;ConfoundTCsNames];
ConfoundTCsNamesRaw=[ConfoundTCsBySubjNames;ConfoundTCsNamesRaw];
if numRuns>1
    ConfoundTCs=[ConfoundTCs,runTC];
    ConfoundTCsNames=[ConfoundTCsNames;runTCNames]; 
    numConfoundTCs=size(ConfoundTCs,2);
else
    numConfoundTCs=size(ConfoundTCs,2)+1;
end    

numEvents=size(Events,2);
numTimeCourses=size(TimeCourses,2);
XTypes=[ones(1,numEvents),ones(1,numTimeCourses)*2,ones(1,numConfoundTCs)*3];
if strcmpi(Normalize,'zscore')
    Events=zscore(Events);
    TimeCourses=zscore(TimeCourses);
    %CounfoundTCs=zscore(ConfoundTCs);
end
X=single([Events,TimeCourses,ConfoundTCs]);
XNames=[EventNames;TimeCoursesNames;ConfoundTCsNames];
XNamesRaw=[EventNames;TimeCoursesNames;ConfoundTCsNamesRaw];
[~,ModeFilter]=mode(X,1);
ModeFilter=ModeFilter==size(X,1);
if numRuns==1
    X=[X,ones(size(X,1),1)];
    XNames=[XNames;{'Constant'}];
    XNamesRaw=[XNamesRaw;{'Constant'}];
    EmptyRegressor=[ModeFilter,0];
else
    EmptyRegressor=ModeFilter;
end

end

function [X] = HRFConvolve(X)
    HRF_use=hrf('twogamma',0.5);
    for i = 1:size(X,2)
        tempTC=conv(X(:,i),HRF_use);
        X(:,i)=tempTC(1:size(X,1),1);
    end
end


function [AlignedVars,AllVarNames,runTC,runTCNames]=NameAlignVar(Vars,VarNames)    
    numRuns = size(VarNames,1);
    AllVarNames=[];
    for run = 1:numRuns
        AllVarNames=unique([AllVarNames;VarNames{run,1}(:)]);
    end
    
    runTC=[];
    AlignedVars=[];
    runTCNames=[];
    for run=1:numRuns  
        runDur=size(Vars{run,1},1);
        tempRunTC=zeros(runDur,numRuns);
        tempRunTC(:,run)=ones(runDur,1);
        runTC=[runTC;tempRunTC];
        runTCNames=[runTCNames;{['Run_',num2str(run)]}];
        [regNameInd,regNameOrder]=ismember(AllVarNames,VarNames{run,1});
        tempVars=zeros(runDur,length(AllVarNames));
        tempVars(:,regNameInd)=Vars{run,1}(:,regNameOrder(regNameInd));
        AlignedVars=[AlignedVars;tempVars];                                      
    end      
end

function [ConfoundMat,ConfoundNames,ConfoundNamesRaw]=cell2ConfoundMat(ConfoundCell,Normalize)
    numRuns=size(ConfoundCell,1);
    matSizes=zeros(numRuns,2);
    ConfoundNames=[];
    ConfoundNamesRaw=[];
    for run = 1:numRuns
        matSizes(run,:)=size(ConfoundCell{run,1});
        tempNames=ConfoundCell{run,2};
        for i = 1:length(tempNames)
            tempNames{i,1}=[tempNames{i,1},'_r',num2str(run)];
        end
        ConfoundNames=[ConfoundNames;tempNames];
        ConfoundNamesRaw=[ConfoundNamesRaw;ConfoundCell{run,2}];
    end
    ConfoundMat=zeros(sum(matSizes,1));    
    endInds=cumsum(matSizes);
    startInds=[1,1;(endInds(1:end-1,:)+1)];
    
    for run = 1:numRuns
        if strcmpi(Normalize,'zscore')
            ConfoundMat(startInds(run,1):endInds(run,1),startInds(run,2):endInds(run,2))=zscore(ConfoundCell{run,1});
        else
            ConfoundMat(startInds(run,1):endInds(run,1),startInds(run,2):endInds(run,2))=ConfoundCell{run,1};
        end
    end
end

