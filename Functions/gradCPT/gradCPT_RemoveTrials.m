function [RemoveTrials,tcInfo] = gradCPT_RemoveTrials( subInfo,BehaviorDataPrefix,Mountains,OEs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

numSS=subInfo.numSS;
numRuns=subInfo.numRuns;
RemoveTrials=cell(1,sum(cell2mat(numRuns),1));
count=1;
for n=1:numSS 
    subID=subInfo.subID{n,1};
    subInfo.runID{n,1}=subInfo.runID{n,1}(:);
    for r = 1:numRuns{n,1}
        runID=subInfo.runID{n,1}{r,1};
        loadName=[BehaviorDataPrefix,subID,runID,subInfo.behavior.suffix];
        [~,~,tcInfo{1,count} ] = gradCPT_TCsByItem(loadName,'full',0);
        RemoveTrials{1,count}=zeros(length(tcInfo{1,count}.OEs{1,1}),1);
        if OEs==1
            RemoveTrials{1,count}=RemoveTrials{1,count}+single(tcInfo{1,count}.OEs{1,1});
        elseif Mountains==1
            RemoveTrials{1,count}=RemoveTrials{1,count}+single(tcInfo{1,count}.StimType{1,1}~=2);            
        end
        count=count+1;
    end
end

end

