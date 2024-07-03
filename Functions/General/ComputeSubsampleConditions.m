function [SubConditions,SubConditionNames,SplitAssigns,SplitData,SplitAssignsByTime] = ComputeSubsampleConditions(Conditions,ConditionNames,SplitVar,CurrentSub,TotalSub,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [SubSampleSize] = VariableSetter('SubSampleSize',0.5,varargin);
    [TrialNum] = VariableSetter('TrialNum',[1:size(Conditions,1)]',varargin);
    SplitVar=SplitVar(:);
    if length(SplitVar)~=size(Conditions,1)
        SplitVar=imresize(SplitVar,[size(Conditions,1),1],'nearest');
    end    
    MaxSub=(TotalSub-1)/2;
    CurrentSub=CurrentSub-(MaxSub+1);
    numCond=length(ConditionNames);
    SubConditions=[Conditions,Conditions]*0;
    SubConditions(SubConditions==0)=nan;
    SplitAssignsByTime=TrialNum;
    NanInd=isnan(TrialNum);
    SplitVar(NanInd,:)=[];
    TrialNum(NanInd,:)=[];
    Conditions(NanInd,:)=[];
    ConditionNames=ConditionNames(:);
    [~,UniqueInd]=unique(TrialNum);
    UniqueConditions=Conditions(UniqueInd,:);
    UniqueSplitVar=SplitVar(UniqueInd,:);
    [ SplitAssigns,SplitData ] = GetSplitDiff(UniqueSplitVar,[CurrentSub,MaxSub],'Conditions',UniqueConditions,'SubSampleSize',SubSampleSize,'reps',1000);
    SplitAssgnsMat=repmat(SplitAssigns,[1,size(UniqueConditions,2)]);
    tempSplitConditions=[UniqueConditions.*single(SplitAssgnsMat==1),UniqueConditions.*single(SplitAssgnsMat==2)];
    tempSplitConditions=tempSplitConditions(TrialNum,:);
%     tempSplitAssignsByTime=SplitAssigns(TrialNum,:);
    SplitAssignsByTime(NanInd==0,1)=SplitAssigns(TrialNum,:);
    SubConditions(NanInd==0,:)=tempSplitConditions;
    SubConditions=[SubConditions(:,1:numCond),sum(SubConditions(:,numCond+1:end),2)];
    SubConditionNames=[ConditionNames;'subRemainder'];
    
    
end

