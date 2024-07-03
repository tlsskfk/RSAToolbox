function [SplitConditions,SplitConditionNames,SplitAssigns,SplitData] = ComputeSplitConditions(Conditions,ConditionNames,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [SplitVar] = VariableSetter('SplitVar',[],varargin);
    [TrialNum] = VariableSetter('TrialNum',[1:size(Conditions,1)]',varargin);
    [RandomSplit] = VariableSetter('RandomSplit',0,varargin);
    if isempty(SplitVar)
        SplitVar=rand(size(Conditions,1),1);
    end
    if length(SplitVar)~=size(Conditions,1)
        SplitVar=imresize(SplitVar,[size(Conditions,1),1],'nearest');
    end
    SplitConditions=[Conditions,Conditions]*0;
    SplitConditions(SplitConditions==0)=nan;
    NanInd=isnan(TrialNum);
    SplitVar(NanInd,:)=[];
    TrialNum(NanInd,:)=[];
    Conditions(NanInd,:)=[];
    ConditionNames=ConditionNames(:);
    [~,UniqueInd]=unique(TrialNum);
    UniqueConditions=Conditions(UniqueInd,:);
    UniqueSplitVar=SplitVar(UniqueInd,:);
    [ tempSplitAssigns,SplitData ] = GetSplitDiff(UniqueSplitVar,[1,1],'Conditions',UniqueConditions,'SubSampleSize',0.5,'reps',10000);
    if RandomSplit==0
        SplitData.groupMeans=flip(SplitData.groupMeans,2);
        SplitData.RealDiff=SplitData.RealDiff*-1;
        SplitData.byCond.groupMeans=flip(SplitData.byCond.groupMeans,2);
        SplitData.byCond.RealDiff=SplitData.byCond.RealDiff*-1;
        SplitAssigns=tempSplitAssigns*0;
        SplitAssigns(tempSplitAssigns==1)=2;
        SplitAssigns(tempSplitAssigns==2)=1;
    else
        SplitAssigns=tempSplitAssigns;
    end
    SplitAssgnsMat=repmat(SplitAssigns,[1,size(UniqueConditions,2)]);
    tempSplitConditions=[UniqueConditions.*single(SplitAssgnsMat==1),UniqueConditions.*single(SplitAssgnsMat==2)];
    tempSplitConditions=tempSplitConditions(TrialNum,:);
    SplitConditions(NanInd==0,:)=tempSplitConditions;
    ConditionNames2=ConditionNames;
    for i = 1:length(ConditionNames)
        ConditionNames{i,1}=[ConditionNames{i,1},'_1'];
        ConditionNames2{i,1}=[ConditionNames2{i,1},'_2'];
    end
    SplitConditionNames=[ConditionNames;ConditionNames2];
end

