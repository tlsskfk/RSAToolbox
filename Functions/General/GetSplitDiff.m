function [ SplitAssigns,SplitData ] = GetSplitDiff(SplitVar,Diff,varargin)
%Written by David Rothlein
%Functions called: VariableSetter, vararginConvert, randSort, Split2GroupsByDiff
%Input: 
%SplitVar is a 1D array of values used to split the data into 2 groups
%Diff is a (1 by 2) array [I,N] where N is an integer indicating the total
%number of splits and I indicates this particular split


SplitVar=SplitVar(:);%Ensure SplitVar is a vertical array
Diff=Diff(:)'; %Ensure Diff is a horizontal array
negDiff=0;
if Diff(1,1) < 0
    SplitVar=SplitVar*-1;
    Diff(1,1)=abs(Diff(1,1));
    negDiff=1;
end
%Set optional variables
%Defaults:
% Conditions = ones(length(SplitVar),1);
% SubSampleSize = 0.5;
% reps = 10000;
% ErrorThresh = 0;

Conditions = VariableSetter( 'Conditions',ones(length(SplitVar),1),varargin);
SubSampleSize = VariableSetter( 'SubSampleSize',0.5,varargin);
reps = VariableSetter( 'reps',10000,varargin);
ErrorThresh = VariableSetter( 'ErrorThresh',0,varargin);

SplitAssigns=zeros(length(SplitVar),1);
numCond=size(Conditions,2);
SplitData.byCond.RealDiff=nan(numCond,1);
SplitData.byCond.Error=nan(numCond,1);
SplitData.byCond.groupMeans=nan(numCond,2);
SplitVar(isinf(SplitVar))=nan;

for cond=1:numCond
    %Define condition specific index, used to pull condition specific splitvar vals 
    %from full timecourse and place condition specific split assigns 
    %back into full timecourse
    condInd=find(Conditions(:,cond)==1);
    Vals=SplitVar(condInd,:);
    
    %Clean by removing NANs and INFs from splitvar vals and also removing
    %those timepoints from the condition index. Such timepoints will not be
    %assigned a split (hence ignored or treated as baseline)
    cleanInd=isnan(Vals);
    Vals(cleanInd,:)=[];
    condInd(cleanInd,:)=[];
    numVals=length(Vals);
      
    
    %It is possible that the splitvar vals do not vary (e.g. all 1 or all 0) in
    %such a case, the split assigns will be allocated randomly proprtional
    %to subsample size.  
    if length(unique(Vals))==1
        SplitAssign=ones(length(Vals),1)+1; 
        rInd=randperm(numVals);        
        %In instances where the length is odd, the majority is allocated to
        %a split at random
        if round(rand(1))==1
            rInd=rInd(1,1:ceil(numVals*SubSampleSize));
        else
            rInd=rInd(1,1:floor(numVals*SubSampleSize));
        end
        SplitAssign(rInd)=1;
        SplitAssigns(condInd)=SplitAssign;
        SplitData.byCond.RealDiff(cond,1)=0;
        SplitData.byCond.groupMeans(cond,1)=[nanmean(Vals(:)),nanmean(Vals(:))];  
        if negDiff == 1
            SplitData.byCond.groupMeans(cond,1)=SplitData.byCond.groupMeans(cond,1)*-1;
        end    
        %skip to next condition iteration
        continue
    end
    
    %Sort splitvar vals, rank of equal values are randomly permuted.
    sortInd=flipud(randSort(Vals)); %Fuction randSort
    maxGroup=zeros(numVals,1);
    if round(rand(1))==1
        maxGroup(sortInd(1:ceil(numVals*SubSampleSize),1),1)=1;
    else
        maxGroup(sortInd(1:floor(numVals*SubSampleSize),1),1)=1;
    end
    maxDiff=nanmean(Vals(maxGroup==1),1)-nanmean(Vals(maxGroup==0),1);
    
    %if Diff var indicates max split, set Group assigns as max diff,
    %else identify the target Diff value and run the gradiant search to
    %split Vals into two groups that differ by target Diff
    if abs(Diff(1,1))==Diff(1,2)
        SplitAssign=ones(length(Vals),1)+1; 
        if negDiff == 1
            SplitData.byCond.RealDiff(cond,1)=maxDiff*-1;
        else    
            SplitData.byCond.RealDiff(cond,1)=maxDiff;
        end
        SplitData.byCond.Error(cond,1)=0;
        SplitAssign(maxGroup==1,1)=1;
        SplitAssigns(condInd)=SplitAssign;
        SplitData.byCond.groupMeans(cond,:)=[nanmean(Vals(maxGroup==1),1),nanmean(Vals(maxGroup==0),1)];
        continue
    elseif Diff(1,1)==0
        TargetDiff=0;
    else
        diffSteps=[maxDiff/(Diff(1,2)):maxDiff/(Diff(1,2)):maxDiff];
        TargetDiff=diffSteps(1,(Diff(1,1)));
    end
    [ SplitAssign,SplitData.byCond.groupMeans(cond,:),SplitData.byCond.RealDiff(cond,1),SplitData.byCond.Error(cond,1) ] = Split2GroupsByDiff( Vals,TargetDiff,SubSampleSize,reps,ErrorThresh );
    if negDiff == 1
        SplitData.byCond.groupMeans(cond,:)=SplitData.byCond.groupMeans(cond,:)*-1;
        SplitData.byCond.RealDiff(cond,1)=SplitData.byCond.RealDiff(cond,1)*-1;
    end
    SplitAssign=single(SplitAssign);
    SplitAssign(SplitAssign==0)=2;
    SplitAssigns(condInd)=SplitAssign;
end   
if negDiff == 1
    SplitVar=SplitVar*-1;
end
SplitData.groupMeans=[nanmean(SplitVar(SplitAssigns==1),1),nanmean(SplitVar(SplitAssigns==2),1)];
SplitData.RealDiff=SplitData.groupMeans(1,2)-SplitData.groupMeans(1,1);
SplitData.Error=nanmean(SplitData.byCond.Error(cond,1),1);
end
