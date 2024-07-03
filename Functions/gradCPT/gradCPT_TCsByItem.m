function [ TCsbyRun,TCbySS,SaveInfo ] = gradCPT_TCsByItem(FilePaths,Type,noOE,trialDur,setDiff) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    trialDur=0.8;
    setDiff=[];
end
if nargin==4
    setDiff=[];
end
SaveInfo.trialDur=trialDur;
SaveInfo.Type=Type;
SaveInfo.noOE=noOE;
if ~iscell(FilePaths)
    FilePaths={FilePaths};
end

numRuns=length(FilePaths);
if size(FilePaths,1)==numRuns
    FilePaths=FilePaths';
end
TCbySS=[];

for i = 1:numRuns
    bData=load(FilePaths{1,i});
    StimType=bData.data(:,3);
    StimType(end,:)=[];
    StimNums=bData.data(:,4);
    StimNums(end,:)=[];
    StimError=bData.response(:,7);
    StimError(end,:)=[];
    OEs=StimType==2 & StimError == 0;
    Mountains=StimType==1;
    SaveInfo.StimType{1,i}=StimType;
    SaveInfo.StimNum{1,i}=StimNums;
    SaveInfo.StimError{1,i}=StimError;
    SaveInfo.OEs{1,i}=OEs;
    SaveInfo.Mountains{1,i}=Mountains;
    SaveInfo.startTime{1,i}=bData.starttime;
    SaveInfo.trialOnsets{1,i}=bData.data(1:end-1,9);
    SaveInfo.localEventOnsets{1,i}=SaveInfo.trialOnsets{1,i}-SaveInfo.startTime{1,i}+(trialDur/2);
    SaveInfo.Duration=bData.endtime-bData.starttime;
    if strcmpi(Type,'zone')==1
        VTC=CPT_analyze_zone_func2(bData.response,bData.data);
        [~,TCsbyRun{1,i},~,~,~] = ZoneByItem( VTC,StimType,StimNums );
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}];   
    elseif strcmpi(Type,'zonePrime')==1
        VTC=CPT_analyze_zone_func2(bData.response,bData.data);
        [~,~,~,TCsbyRun{1,i},~] = ZoneByItem( VTC,StimType,StimNums );
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'reward')==1
        TCsbyRun{1,i}=bData.bordertracker(:,2)==255;
        TCsbyRun{1,i}(end,:)=[];
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'full')
        TCsbyRun{1,i}=ones(length(StimType),1);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}];   
    elseif strcmpi(Type,'zoneSample')==1
        VTC=CPT_analyze_zone_func2(bData.response,bData.data);
        [~,TCsbyRun{1,i},~,~,~] = ZoneByItem( VTC,StimType,StimNums,setDiff );
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'zoneDerivSample')==1
        VTC=CPT_analyze_zone_func2(bData.response,bData.data);
        [~,~,~,TCsbyRun{1,i},~] = ZoneByItem( VTC,StimType,StimNums,setDiff );
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'rewardSample')==1
        TCsbyRun{1,i}=bData.bordertracker(:,2)==255;
        TCsbyRun{1,i}(end,:)=[];
        [ TCsbyRun{1,i},~,~ ] = FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    else
        SaveInfo.Type='Rand';
        VTC=CPT_analyze_zone_func2(bData.response,bData.data);
        VTC=VTC(randperm(length(VTC)),:);
        [~,TCsbyRun{1,i},~,~,~] = ZoneByItem( VTC,StimType,StimNums );
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    end
        
end

