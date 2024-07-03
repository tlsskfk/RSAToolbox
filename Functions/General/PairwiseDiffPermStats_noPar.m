function [Results] = PairwiseDiffPermStats_noPar(InData,NumReps,fixedVal,tails)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    NumReps=10000;
    fixedVal=0;
    tails='two-tailed';
end
if nargin == 2
    fixedVal=0;
    tails='two-tailed';
end
if nargin == 3
    tails='two-tailed';
end

numSets=size(InData,2);
numGroups=size(InData,1);
Results.PairwisePermPs=nan(numGroups,numGroups,numSets);
Results.PairwiseTPs=nan(numGroups,numGroups,numSets);
Results.PairwiseRealDiff=nan(numGroups,numGroups,numSets);
Results.PairwiseTVal=nan(numGroups,numGroups,numSets);
Results.PairwisePermCI=cell(numGroups,numGroups,numSets);
Results.PairwiseTCI=cell(numGroups,numGroups,numSets);

IndvPermPs=nan(numGroups*numSets,1);
IndvTPs=nan(numGroups*numSets,1);
IndvRealDiff=nan(numGroups*numSets,1);
IndvTVal=nan(numGroups*numSets,1);
IndvPermCI=cell(numGroups*numSets,1);
IndvTCI=cell(numGroups*numSets,1);
tempInData=InData(:);
for i = 1:length(tempInData)
    [IndvPermPs(i,1),IndvTPs(i,1),IndvRealDiff(i,1),IndvPermCI{i,1},IndvTCI{i,1},IndvTVal(i,1)] = PermDiffFcn(tempInData{i,1},[],NumReps,fixedVal,tails);
end
Results.IndvPermPs=reshape(IndvPermPs,[numGroups,numSets]);
Results.IndvTPs=reshape(IndvTPs,[numGroups,numSets]);
Results.IndvRealDiff=reshape(IndvRealDiff,[numGroups,numSets]);
Results.IndvTVal=reshape(IndvTVal,[numGroups,numSets]);
Results.IndvPermCI=reshape(IndvPermCI,[numGroups,numSets]);
Results.IndvTCI=reshape(IndvTCI,[numGroups,numSets]);

if numGroups>1
    AllPairs=[];
    for i = 1:numSets
        tempPairs=nchoosek([1:numGroups],2);
        AllPairs=[AllPairs;[tempPairs,ones(size(tempPairs,1),1)*i]];
    end
    numPairs=size(AllPairs,1);
    tempPairwisePermPs=cell(numPairs,1);
    tempPairwiseTPs=cell(numPairs,1);
    tempPairwiseRealDiff=cell(numPairs,1);
    tempPairwisePermCI=cell(numPairs,1);
    tempPairwiseTCI=cell(numPairs,1);
    tempPairwiseTVal=cell(numPairs,1);    
    for j = 1:size(AllPairs,1)
        Data1Ind=AllPairs(j,1);
        Data2Ind=AllPairs(j,2);
        setInd=AllPairs(j,3);
        [tempPairwisePermPs{j,1},tempPairwiseTPs{j,1},tempPairwiseRealDiff{j,1},tempPairwisePermCI{j,1},tempPairwiseTCI{j,1},tempPairwiseTVal{j,1}] = PermDiffFcn(InData{Data1Ind,setInd},InData{Data2Ind,setInd},NumReps,fixedVal,tails);
    end
    for j = 1:size(AllPairs,1)
        Data1Ind=AllPairs(j,1);
        Data2Ind=AllPairs(j,2);
        setInd=AllPairs(j,3);
        Results.PairwisePermPs(Data1Ind,Data2Ind,setInd)=tempPairwisePermPs{j,1};
        Results.PairwiseTPs(Data1Ind,Data2Ind,setInd)=tempPairwiseTPs{j,1};
        Results.PairwiseRealDiff(Data1Ind,Data2Ind,setInd)=tempPairwiseRealDiff{j,1};
        Results.PairwisePermCI{Data1Ind,Data2Ind,setInd}=tempPairwisePermCI{j,1};
        Results.PairwiseTCI{Data1Ind,Data2Ind,setInd}=tempPairwiseTCI{j,1};
        Results.PairwiseTVal(Data1Ind,Data2Ind,setInd)=tempPairwiseTVal{j,1};        
    end
else
    PairwisePermPs=[];
    PairwiseTPs=[];
end
        
       
end

function [PermP,tP,realDiff,permCI,tCI,tVal] = PermDiffFcn(Data1,Data2,NumReps,fixedVal,tails)
    bsDist=zeros(NumReps,1);
    if isempty(Data2)
        for bsRep=1:NumReps
            bsDist(bsRep,1)=mean(Data1(randi(length(Data1),[length(Data1),1]),1));
        end
        [tempP] = reverseprctile(bsDist,fixedVal);
        if strcmpi(tails,'two-tailed')
            if tempP > 0.5
                tempP=1-tempP;
            end
            PermP=tempP*2;
        elseif strcmpi(tails,'left-tailed')
            PermP=1-tempP;
        else
            PermP=tempP;
        end
        tCI=zeros(1,2);
        permCI=prctile(bsDist,[2.5,97.5]);
        [tVal,~,realDiff,tCI(1,1),tCI(1,2),tP]=getTval(Data1,1);
    else
        realDiff=nanmean(Data1,1)-nanmean(Data2,1);
        PullInd=size(Data1,1);
        AllData=[Data1;Data2];
        for bsRep=1:NumReps
            tempAllData=AllData(randperm(length(AllData)),:);
            bsDist(bsRep,1)=nanmean(tempAllData(1:PullInd,1),1)-nanmean(tempAllData(PullInd+1:end,1),1);
        end
        
        if strcmpi(tails,'two-tailed')
            PermP = 1 - reverseprctile(abs(bsDist),abs(realDiff));
            permCI=prctile(bsDist,[2.5,97.5]);
            [~,tP,tCI,tempStats]=ttest2(Data1,Data2);
            tVal=tempStats.tstat;
        elseif strcmpi(tails,'left-tailed')
            PermP=reverseprctile(bsDist,realDiff);
            permCI=prctile(bsDist,[2.5,97.5]);
            [~,tP,tCI,tempStats]=ttest2(Data1,Data2,'tail','left');
            tVal=tempStats.tstat;
        else
            tempP=reverseprctile(bsDist,realDiff);
            PermP=1-tempP;
            permCI=prctile(bsDist,[2.5,97.5]);
            [~,tP,tCI,tempStats]=ttest2(Data1,Data2,'tail','right');
            tVal=tempStats.tstat;
        end
    end  
    tCI=tCI';
end


