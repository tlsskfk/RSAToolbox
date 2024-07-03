function [ TCsbyRun,TCbySS,SaveInfo ] = gradCPT_TCsByItem2(TrialData,Type,noOE,trialDur,setDiff) 
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

numRuns=length(TrialData);

TCbySS=[];

for i = 1:numRuns
    SaveInfo.StimType{1,i}=TrialData{i,1}.StimTypes;
    SaveInfo.StimNum{1,i}=TrialData{i,1}.StimNums;
    SaveInfo.StimError{1,i}=TrialData{i,1}.Events.CE;
    SaveInfo.OEs{1,i}=TrialData{i,1}.Events.OE;
    SaveInfo.Mountains{1,i}=TrialData{i,1}.Events.mountain;
    %SaveInfo.startTime{1,i}=bData.starttime;
    %SaveInfo.trialOnsets{1,i}=bData.data(1:end-1,9);
    SaveInfo.localEventOnsets{1,i}=TrialData{i,1}.TrialOnset;
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
        TCsbyRun{1,i}=TrialData{i,1}.TimeCourses.VTC;
        %[Zone,ZonebyItem,ZonePrime,ZonePrimebyItem,VTCderiv] = ZoneByItem( VTC,StimType,StimNums,setDiff )
        [~,TCsbyRun{1,i},~,~,~] = ZoneByItem( TCsbyRun{1,i},TrialData{i,1}.StimTypes,TrialData{i,1}.StimNums,setDiff );
        %TCsbyRun{1,i}=FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'zoneDerivSample')==1
        TCsbyRun{1,i}=TrialData{i,1}.TimeCourses.VTCderiv;
        TCsbyRun{1,i}=FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'rewardSample')==1
        TCsbyRun{1,i}=TrialData{i,1}.TimeCourses.Reward;
        TCsbyRun{1,i}=FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'mwSample')==1
        TCsbyRun{1,i}=TrialData{i,1}.TimeCourses.MindWanderTP;
        TCsbyRun{1,i}=FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
        if noOE==1
            TCsbyRun{1,i}(OEs,1)=-1;
        end
        TCbySS=[TCbySS;TCsbyRun{1,i}]; 
    elseif strcmpi(Type,'confidenceSample')==1
        TCsbyRun{1,i}=TrialData{i,1}.TimeCourses.ConfidenceTP;
        TCsbyRun{1,i}=FindGroupDiff(single(TCsbyRun{1,i}),setDiff,10000);
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

