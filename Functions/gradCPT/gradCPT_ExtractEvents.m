function [ fmriprep_table ] = gradCPT_ExtractEvents(fmriprep_table,ExperimentsDir,varargin)
%EXAMPLE
% filepaths={'/Volumes/BALL_LAB/gradCPT_MindWandering/DataFromEve/MW001_HC/fMRI_Behavior/MW001_HC_Behavior_gradCPT_TP1.mat','/Volumes/BALL_LAB/gradCPT_MindWandering/DataFromEve/MW001_HC/fMRI_Behavior/MW001_HC_Behavior_gradCPT_TP2.mat','/Volumes/BALL_LAB/gradCPT_MindWandering/DataFromEve/MW001_HC/fMRI_Behavior/MW001_HC_Behavior_gradCPT_TP3.mat','/Volumes/BALL_LAB/gradCPT_MindWandering/DataFromEve/MW001_HC/fMRI_Behavior/MW001_HC_Behavior_gradCPT_TP4.mat'};
%[ TrialData ] = gradCPT_ExtractEvents( filepaths,'expType','tp');

stimTypes = VariableSetter( 'stimTypes','citymountain10',varargin);
respTypes = VariableSetter( 'respTypes',[],varargin);
if strcmpi(stimTypes,'citymountain10')
    condNames = VariableSetter( 'condNames',{{0,'TP_Block'},{1,'mountain'},{2,'city'}},varargin);
   % respTypes = VariableSetter( 'respTypes',{{0,'noPress'},{30,'Press'},{21,'Press'}},varargin);
    EventTypes = VariableSetter( 'EventTypes',{{'mountain','Press','CE'},{'mountain','noPress','CO'},{'city','Press','CC'},{'city','noPress','OE'},{'TP_Block','Press','TP_Block'},{'TP_Block','noPress','TP_Block'}},varargin);
    ExemplarTypes = VariableSetter( 'ExemplarTypes',{{'city',1,'city01'},{'city',2,'city02'},{'city',3,'city03'},{'city',4,'city04'},{'city',5,'city05'},{'city',6,'city06'},{'city',7,'city07'},{'city',8,'city08'},{'city',9,'city09'},{'city',10,'city10'},{'mountain',1,'mountain01'},{'mountain',2,'mountain02'},{'mountain',3,'mountain03'},{'mountain',4,'mountain04'},{'mountain',5,'mountain05'},{'mountain',6,'mountain06'},{'mountain',7,'mountain07'},{'mountain',8,'mountain08'},{'mountain',9,'mountain09'},{'mountain',10,'mountain10'}},varargin);
elseif strcmpi(stimTypes,'emodig')
    condNames = VariableSetter( 'condNames',{{0,'TP_Block'},{1,'three'},{2,'nthree'}},varargin);
   % respTypes = VariableSetter( 'respTypes',{{0,'noPress'},{30,'Press'},{21,'Press'}},varargin);
    EventTypes = VariableSetter( 'EventTypes',{{'three','Press','CE'},{'three','noPress','CO'},{'nthree','Press','CC'},{'nthree','noPress','OE'},{'TP_Block','Press','TP_Block'},{'TP_Block','noPress','TP_Block'}},varargin);
    ExemplarTypes = VariableSetter( 'ExemplarTypes',{{'nthree',1,'large1'},{'nthree',2,'large1_h'},{'nthree',3,'large1_p'},{'nthree',4,'large2'},{'nthree',5,'large2_h'},{'nthree',6,'large2_p'},{'nthree',7,'large4'},{'nthree',8,'large4_h'},{'nthree',9,'large4_p'},{'nthree',10,'large5'},{'nthree',11,'large5_h'},{'nthree',12,'large5_p'},{'nthree',13,'large6'},{'nthree',14,'large6_h'},{'nthree',15,'large6_p'},{'nthree',16,'large7'},{'nthree',17,'large7_h'},{'nthree',18,'large7_p'},{'nthree',19,'large8'},{'nthree',20,'large8_h'},{'nthree',21,'large8_p'},{'nthree',22,'large9'},{'nthree',23,'large9_h'},{'nthree',24,'large9_p'},{'nthree',25,'med1'},{'nthree',26,'med1_h'},{'nthree',27,'med1_p'},{'nthree',28,'med2'},{'nthree',29,'med2_h'},{'nthree',30,'med2_p'},{'nthree',31,'med4'},{'nthree',32,'med4_h'},{'nthree',33,'med4_p'},{'nthree',34,'med5'},{'nthree',35,'med5_h'},{'nthree',36,'med5_p'},{'nthree',37,'med6'},{'nthree',38,'med6_h'},{'nthree',39,'med6_p'},{'nthree',40,'med7'},{'nthree',41,'med7_h'},{'nthree',42,'med7_p'},{'nthree',43,'med8'},{'nthree',44,'med8_h'},{'nthree',45,'med8_p'},{'nthree',46,'med9'},{'nthree',47,'med9_h'},{'nthree',48,'med9_p'},{'nthree',49,'small1'},{'nthree',50,'small1_h'},{'nthree',51,'small1_p'},{'nthree',52,'small2'},{'nthree',53,'small2_h'},{'nthree',54,'small2_p'},{'nthree',55,'small4'},{'nthree',56,'small4_h'},{'nthree',57,'small4_p'},{'nthree',58,'small5'},{'nthree',59,'small5_h'},{'nthree',60,'small5_p'},{'nthree',61,'small6'},{'nthree',62,'small6_h'},{'nthree',63,'small6_p'},{'nthree',64,'small7'},{'nthree',65,'small7_h'},{'nthree',66,'small7_p'},{'nthree',67,'small8'},{'nthree',68,'small8_h'},{'nthree',69,'small8_p'},{'nthree',70,'small9'},{'nthree',71,'small9_h'},{'nthree',72,'small9_p'},{'three',1,'large3'},{'three',2,'large3_h'},{'three',3,'large3_p'},{'three',4,'med3'},{'three',5,'med3_h'},{'three',6,'med3_p'},{'three',7,'small3'},{'three',8,'small3_h'},{'three',9,'small3_p'},{'nthree',1,'one'},{'nthree',2,'one'},{'nthree',3,'one'},{'nthree',4,'two'},{'nthree',5,'two'},{'nthree',6,'two'},{'nthree',7,'four'},{'nthree',8,'four'},{'nthree',9,'four'},{'nthree',10,'five'},{'nthree',11,'five'},{'nthree',12,'five'},{'nthree',13,'six'},{'nthree',14,'six'},{'nthree',15,'six'},{'nthree',16,'seven'},{'nthree',17,'seven'},{'nthree',18,'seven'},{'nthree',19,'eight'},{'nthree',20,'eight'},{'nthree',21,'eight'},{'nthree',22,'nine'},{'nthree',23,'nine'},{'nthree',24,'nine'},{'nthree',25,'one'},{'nthree',26,'one'},{'nthree',27,'one'},{'nthree',28,'two'},{'nthree',29,'two'},{'nthree',30,'two'},{'nthree',31,'four'},{'nthree',32,'four'},{'nthree',33,'four'},{'nthree',34,'five'},{'nthree',35,'five'},{'nthree',36,'five'},{'nthree',37,'six'},{'nthree',38,'six'},{'nthree',39,'six'},{'nthree',40,'seven'},{'nthree',41,'seven'},{'nthree',42,'seven'},{'nthree',43,'eight'},{'nthree',44,'eight'},{'nthree',45,'eight'},{'nthree',46,'nine'},{'nthree',47,'nine'},{'nthree',48,'nine'},{'nthree',49,'one'},{'nthree',50,'one'},{'nthree',51,'one'},{'nthree',52,'two'},{'nthree',53,'two'},{'nthree',54,'two'},{'nthree',55,'four'},{'nthree',56,'four'},{'nthree',57,'four'},{'nthree',58,'five'},{'nthree',59,'five'},{'nthree',60,'five'},{'nthree',61,'six'},{'nthree',62,'six'},{'nthree',63,'six'},{'nthree',64,'seven'},{'nthree',65,'seven'},{'nthree',66,'seven'},{'nthree',67,'eight'},{'nthree',68,'eight'},{'nthree',69,'eight'},{'nthree',70,'nine'},{'nthree',71,'nine'},{'nthree',72,'nine'}},varargin);
end    

AddTimeConstant = VariableSetter( 'AddTimeConstant',[0],varargin);
BlockThresh = VariableSetter( 'BlockThresh',8,varargin);
gauss_sigma = VariableSetter( 'gauss_sigma',20,varargin);
expType = VariableSetter( 'expType','',varargin);
ScanDur = VariableSetter( 'ScanDur',[],varargin);
if isempty(respTypes)
    respTypes={{0,'noPress'}};
    for i = 1:200
        respTypes=[respTypes,{{i,'Press'}}];
    end
end
shiftStarttime = VariableSetter( 'shiftStarttime',0,varargin);
Overwrite = VariableSetter( 'Overwrite',[0],varargin);

condNames=condNames(:)';
respTypes=respTypes(:)';
EventTypes=EventTypes(:)';
ExemplarTypes=ExemplarTypes(:)';

numConds=length(condNames);
numEvents=length(EventTypes);
numResp=length(respTypes);
numExemplars=length(ExemplarTypes);
OutVars=struct;
if length(unique(fmriprep_table.session))>1
    useIndicies=find(fmriprep_table.numRuns_bySes)';
    numRuns=fmriprep_table.numRuns_bySes;        
else
    useIndicies=find(fmriprep_table.numRuns_bySub)';
    numRuns=fmriprep_table.numRuns_bySub;
end
for ssInd=useIndicies
    for run=1:numRuns(ssInd,1)
        dataInd=ssInd+run-1;
        disp(dataInd)
        expType=fmriprep_table.task{dataInd,1};
        SaveDir=strrep(fmriprep_table.funcDir{dataInd,1},'/fmriprep/','/matlab/');
        SaveDirClean=strrep(SaveDir,'/func/',['/beh/clean/']);
        SaveDirEvents=strrep(SaveDir,'/func/',['/beh/events/']);
        LoadDir=strrep(SaveDir,'/func/',['/beh/raw/']);
        SaveNameClean=['sub-',fmriprep_table.sub{dataInd,1},'_task-',fmriprep_table.task{dataInd,1},'_run-',num2str(fmriprep_table.run(dataInd,1)),'_desc-beh_clean.mat'];
        SaveNameEvents=['sub-',fmriprep_table.sub{dataInd,1},'_task-',fmriprep_table.task{dataInd,1},'_run-',num2str(fmriprep_table.run(dataInd,1)),'_desc-beh_events.mat'];
        LoadName=['sub-',fmriprep_table.sub{dataInd,1},'_task-',fmriprep_table.task{dataInd,1},'_run-',num2str(fmriprep_table.run(dataInd,1)),'_desc-beh_raw.mat'];
        try
            load([ExperimentsDir,LoadDir,LoadName],'beh_raw');
        catch
            disp(['No data! Skipping- ',ExperimentsDir,LoadDir,LoadName])
            continue
        end
        data=beh_raw.data;
        try
            endtime=beh_raw.endtime;
        catch
            endtime=beh_raw.ttt(end,end);
        end
        try
            starttime=beh_raw.starttime;
            if shiftStarttime ~= 0
                starttime=starttime+shiftStarttime;
            end
        catch
            starttime=beh_raw.ttt(1,1);
        end
        response=beh_raw.response;
        ttt=beh_raw.ttt;        
        beh_clean=struct;
        beh_clean.Duration=diff(data(:,9));

        beh_clean.ExperimentDuration=endtime-starttime;
        beh_clean.StartTimePad=ttt(1,1)-starttime;
        beh_clean.EndTimePad=endtime-ttt(end,1); 
        response2=response;
        response(end,:)=[];
        data2=data;
        data(end,:)=[];
        numTrials=length(data);
        beh_clean.TrialNum=[1:numTrials]';
        %Create matrix to convert trial space into time space at a 100ms
        %resolution
        beh_clean.Trial2Time=zeros(round(beh_clean.ExperimentDuration*10),1);
        try
        beh_clean.Block=single(beh_clean.Duration>=BlockThresh);
        beh_clean.Event=single(beh_clean.Duration<BlockThresh);
        beh_clean.AverageEventDur=nanmean(beh_clean.Duration(beh_clean.Event==1,:),1);
        fmriprep_table.StartTimePad{dataInd,1}=beh_clean.StartTimePad;
        fmriprep_table.EndTimePad{dataInd,1}=beh_clean.EndTimePad;
        fmriprep_table.ExperimentDuration{dataInd,1}=beh_clean.ExperimentDuration;
        fmriprep_table.AverageEventDur{dataInd,1}=beh_clean.AverageEventDur;        
        if exist([ExperimentsDir,SaveDirClean,SaveNameClean],'file') && Overwrite==0            
            continue
        elseif ~exist([ExperimentsDir,SaveDirClean])
            mkdir([ExperimentsDir,SaveDirClean]);
        end
        beh_clean.StimTypes=data(:,3);
        beh_clean.StimNums=data(:,4);
        %Block onset at beginning of trial; Event onset mid trial (at peak
        %stimulus clarity)
        %beh_clean.TrialOnset=data(:,9)-starttime+(beh_clean.Duration.*beh_clean.Event)+AddTimeConstant;
        beh_clean.TrialOnset=data(:,9)-starttime+AddTimeConstant;
        Trial2TimeIndex=[round((ttt(1:end-1,1)-starttime)*10)+1,round((ttt(2:end,1)-starttime)*10)];
        if sum(beh_clean.Block)>0
            VTC=nan(length(beh_clean.Block),1);
            VTCseg=[[1;find(beh_clean.Block)+1],[find(beh_clean.Block)-1;length(beh_clean.Block)]];
            for BlockSeg=1:size(VTCseg,1)
                bLength=VTCseg(BlockSeg,2)-VTCseg(BlockSeg,1);
                maxL=floor((bLength+2)/3);
                if maxL>20
                    maxL=20;
                end
                VTC(VTCseg(BlockSeg,1):VTCseg(BlockSeg,2),1)=CPT_analyze_zone_func2(response2(VTCseg(BlockSeg,1):VTCseg(BlockSeg,2)+1,:),data2(VTCseg(BlockSeg,1):VTCseg(BlockSeg,2)+1,:),maxL);
            end    
        else
            VTC=CPT_analyze_zone_func2(response2,data2);
        end

        for trial=1:numTrials
            beh_clean.Trial2Time(Trial2TimeIndex(trial,1):Trial2TimeIndex(trial,2),1)=trial;
            condNum=data(trial,3);
            respNum=response(trial,2);       
            for condi=1:numConds
                if condNames{1,condi}{1,1}==condNum
                    beh_clean.Condition{trial,1} = condNames{1,condi}{1,2};
                    break
                end
            end
            beh_clean.Response{trial,1} = 'NA';
            for respi=1:numResp
                if respTypes{1,respi}{1,1}==respNum
                    beh_clean.Response{trial,1} = respTypes{1,respi}{1,2};
                    break
                end
            end  

        end
        for condi=1:numConds
            beh_clean.Events.(condNames{1,condi}{1,2})=single(strcmpi(beh_clean.Condition,condNames{1,condi}{1,2}));
        end
        for respi=1:numResp
            beh_clean.Events.(respTypes{1,respi}{1,2})=single(strcmpi(beh_clean.Response,respTypes{1,respi}{1,2}));
        end    

        for eventi=1:numEvents
            temp=single((beh_clean.Events.(EventTypes{1,eventi}{1,1})+beh_clean.Events.(EventTypes{1,eventi}{1,2}))==2);
            if isfield(beh_clean.Events,EventTypes{1,eventi}{1,3})
                beh_clean.Events.(EventTypes{1,eventi}{1,3})=single((beh_clean.Events.(EventTypes{1,eventi}{1,3})+temp)~=0);
            else
                beh_clean.Events.(EventTypes{1,eventi}{1,3})=temp;
            end
        end

        for exempi=1:numExemplars
            temp=single((beh_clean.Events.(ExemplarTypes{1,exempi}{1,1})+single(data(:,4)==ExemplarTypes{1,exempi}{1,2}))==2);
            if exempi~=1
                if isfield(beh_clean.Exemplars,(ExemplarTypes{1,exempi}{1,3}))
                    beh_clean.Exemplars.(ExemplarTypes{1,exempi}{1,3})=single((beh_clean.Exemplars.(ExemplarTypes{1,exempi}{1,3})+temp)~=0);
                else
                    beh_clean.Exemplars.(ExemplarTypes{1,exempi}{1,3})=temp;
                end
            else
                beh_clean.Exemplars.(ExemplarTypes{1,exempi}{1,3})=temp;   
            end
        end  

        beh_clean.TimeCourses.VTC=VTC;
        beh_clean.TimeCourses.VTCderiv=diff([mean(VTC(:));VTC]);
        RTs=response(:,5);
        RTs(RTs==0)=nan;
        RTs=fillmissing(RTs,'linear');
        SlowRTs=RTs;
        SlowRTs(RTs<=median(RTs))=median(RTs);
        SlowRTs=abs(SlowRTs-median(RTs));
        FastRTs=RTs;
        FastRTs(RTs>median(RTs))=median(RTs);
        FastRTs=abs(FastRTs-median(RTs)); 
        beh_clean.TimeCourses.SlowRTs=SlowRTs;
        beh_clean.TimeCourses.FastRTs=FastRTs;
        beh_clean.TimeCourses.RTs=RTs;
        beh_clean.TimeCourses.FastRTsSmooth=smooth(FastRTs,0.1,'loess');
        beh_clean.TimeCourses.SlowRTsSmooth=smooth(SlowRTs,0.1,'loess');
        beh_clean.TimeCourses.RTsSmooth=smooth(RTs,0.1,'loess');
        beh_clean.TimeCourses.RTs_AR1=AR1_Clean(RTs);
        
        

        if ~isempty(expType)
            if strcmpi(expType,'reward') || strcmpi(expType,'gradCPTwithReward')
                beh_clean.TimeCourses.Reward=single(beh_raw.bordertracker(:,2)==255);
                beh_clean.TimeCourses.Reward(end,:)=[];
            end
            if strcmpi(expType,'thoughtprobe') || strcmpi(expType,'gradCPTMW') || strcmpi(expType,'TP')
                tpInd=data2(:,3)==0;
                temp1=nan(length(tpInd),1);
                temp2=nan(length(tpInd),1);
                temp1(tpInd)=beh_raw.TP_Results(:,1);
                temp2(tpInd)=beh_raw.TP_Results(:,4);
                temp1=fillmissing(temp1,'linear','EndValues','nearest');
                temp1(end,:)=[];
                temp2=fillmissing(temp2,'linear','EndValues','nearest');
                temp2(end,:)=[];
                beh_clean.TimeCourses.MindWanderTP=temp1;
                beh_clean.TimeCourses.ConfidenceTP=temp2;
            end      
        end 
        save([ExperimentsDir,SaveDirClean,SaveNameClean],'beh_clean');
        AllNames=[];
        TrialMat=[];
        if exist([ExperimentsDir,SaveDirEvents,SaveNameEvents],'file') && Overwrite==0            
            continue
        elseif ~exist([ExperimentsDir,SaveDirEvents])
            mkdir([ExperimentsDir,SaveDirEvents]);
        end        
        try
            if iscell(fmriprep_table.TR)
                ScanDur=fmriprep_table.TR{dataInd,1}*fmriprep_table.numVol{dataInd,1};
            else    
                ScanDur=fmriprep_table.TR(dataInd,1)*fmriprep_table.numVol{dataInd,1};
            end
            TimeIncrementInSec=0.5;
        catch
            disp('Cant compute Scan Duration')
            continue
        end
        gauss_sigma=round(20/beh_clean.AverageEventDur); %~20seconds
        ShiftArray=[{-4},{-3},{-2},{-1},{1},{2},{3},{4};{'_L4'},{'_L3'},{'_L2'},{'_L1'},{'_R1'},{'_R2'},{'_R3'},{'_R4'}];
        if isfield(beh_raw,'picNameArray') 
            if ~isfield(beh_raw,'backtime') 
                [beh_raw.backtime] = ttt2backtime(ttt);
            end
            picNameArray=strrep(beh_raw.picNameArray,'.jpg','');
            bg2trialVec=bg2trial(ttt,beh_raw.backtime); 
            bg2trialVec(end,:)=[];
            picNameArray=picNameArray(1:size(bg2trialVec,2),:);
            posArray=contains(picNameArray,'pos');
            negArray=contains(picNameArray,'neg');
            neutArray=~posArray & ~negArray;
            beh_clean.Events.posBG=single(sum(bg2trialVec(:,posArray),2)~=0);    
            for k=1:size(ShiftArray,2)
                VarName=['posBG',ShiftArray{2,k}];
                [beh_clean.Events.(VarName)] = TCShift(beh_clean.Events.posBG,ShiftArray{1,k});
            end    
            beh_clean.Events.negBG=single(sum(bg2trialVec(:,negArray),2)~=0);
            for k=1:size(ShiftArray,2)
                VarName=['negBG',ShiftArray{2,k}];
                [beh_clean.Events.(VarName)] = TCShift(beh_clean.Events.negBG,ShiftArray{1,k});
            end    
            beh_clean.Events.neutBG=single(sum(bg2trialVec(:,neutArray),2)~=0); 
            for k=1:size(ShiftArray,2)
                VarName=['neutBG',ShiftArray{2,k}];
                [beh_clean.Events.(VarName)] = TCShift(beh_clean.Events.neutBG,ShiftArray{1,k});
            end
            bgUniqueNames=unique(picNameArray);
            beh_clean.bgExemplarNames=cell(size(bgUniqueNames,1),1);
            for i = 1:size(bgUniqueNames,1)
                beh_clean.bgExemplarNames{i,1}=['bg_',char(bgUniqueNames(i,1))];
                beh_clean.Exemplars.(beh_clean.bgExemplarNames{i,1})=single(sum(bg2trialVec(:,ismember(picNameArray,bgUniqueNames(i,1))),2)~=0);       
            end
            beh_clean.bgPicNameArray=cell(size(picNameArray,1),1);
            for i = 1:size(picNameArray,1)
                beh_clean.bgPicNameArray{i,1}=char(picNameArray(i,1));
            end
        end    

        if isfield(beh_clean,'Events')
            AllEventNames=fieldnames(beh_clean.Events);
            for NameNum = 1:length(AllEventNames)
                eName=AllEventNames{NameNum,1};
                [ beh_clean.Events.([eName,'_Density']),beh_clean.Events.([eName,'_Pre']),beh_clean.Events.([eName,'_Post']) ] = EventDensityTCs( beh_clean.Events.(eName),0,gauss_sigma);
            end                
        end
        if isfield(beh_clean,'Exemplars')
            AllEventNames=fieldnames(beh_clean.Exemplars);
            for NameNum = 1:length(AllEventNames)
                eName=AllEventNames{NameNum,1};
                [ beh_clean.Exemplars.([eName,'_Density']),beh_clean.Exemplars.([eName,'_Pre']),beh_clean.Exemplars.([eName,'_Post']) ] = EventDensityTCs( beh_clean.Exemplars.(eName),0,gauss_sigma);
            end                
        end        
        
        if isfield(beh_clean,'Events')
            AllNames=[AllNames;fieldnames(beh_clean.Events)];
            TrialMat=[TrialMat,struct2array(beh_clean.Events)];
        end
        if isfield(beh_clean,'Exemplars')
            AllNames=[AllNames;fieldnames(beh_clean.Exemplars)];
            TrialMat=[TrialMat,struct2array(beh_clean.Exemplars)];
        end 
        if isfield(beh_clean,'TimeCourses')
            AllNames=[AllNames;fieldnames(beh_clean.TimeCourses)];
            TrialMat=[TrialMat,struct2array(beh_clean.TimeCourses)];
        end     
        if ~isempty(AllNames)
            AllNames=[AllNames;{'TrialNum'};'TrialOnsetTime'];
            event_names=AllNames;
            TrialMat=[TrialMat,beh_clean.TrialNum,beh_clean.TrialOnset];        
            [resampledTrialMat]=trial2time(TrialMat,beh_clean.Trial2Time,ScanDur,TimeIncrementInSec);
            beh_events=table;
            for i = 1:size(resampledTrialMat,2)
                beh_events=[beh_events,table(resampledTrialMat(:,i),'VariableNames',AllNames(i,1))];
            end
            save([ExperimentsDir,SaveDirEvents,SaveNameEvents],'beh_events');
        end        
        catch
            disp(['Data Error: ',LoadDir,LoadName])
        end
    end
end 

end

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

function bg2trialVec=bg2trial(trialTimes,bgTimes)
    bgFrames=sum(single(bgTimes~=0),2);
    trialFrames=sum(single(trialTimes~=0),2);
    bgFrames(bgFrames==0,:)=[];
    trialFrameNum=unique(trialFrames);
    if length(trialFrameNum)>1
        disp('Warning! Inconsistent trial frame #. Check ttt var');
        trialFrameNum=trialFrameNum(1,1);
        disp(['Using value: ', num2str(trialFrameNum)]);
    end        
    trialVecLengths=bgFrames(:)/trialFrameNum;
    bg2trialVec=zeros(sum(trialVecLengths,1),length(trialVecLengths));
    %count=1;
    count=1+2;
	for i = 1:length(trialVecLengths)
        %vecLength=trialVecLengths(i,1)*2;
        vecLength=trialVecLengths(i,1);
        if count+vecLength-1 >= size(bg2trialVec,1)
            bg2trialVec(count:end,i)=1;
        else
            bg2trialVec(count:count+vecLength-1,i)=1;
        end
        %count=count+vecLength/2;
        count=count+vecLength;
    end    
end