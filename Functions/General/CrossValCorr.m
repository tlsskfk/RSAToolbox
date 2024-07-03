function [ rVals, stats, rperm ] = CrossValCorr( Data1,Data2,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%  'CrossValType': ['leave1out','splithalf']
%  'filterVal': [X]
%  'Reps': [X]
%  'corrType': ['pearson','spearman','kendall']
%  'Norm': ['zScore','fisherZ','both']
%  'whichstats': {'rmean','fisherZmean','tstat','prctiles'}
%  'Permlist': [permMat]  
%  'Tail': ['right' < m>0 > ,'left' < m<0 >,'both' < m~=0 >]

Type = VariableSetter( 'CrossValType','leave1out',varargin);
filterVal = VariableSetter( 'filterVal',[],varargin);
Reps = VariableSetter( 'Reps',1000,varargin);
corrType = VariableSetter( 'corrType','pearson',varargin);
whichstats = VariableSetter( 'whichstats',{''},varargin);
Norm = VariableSetter( 'Norm','',varargin);
rperm = VariableSetter( 'Permlist',[],varargin);
Tail = VariableSetter( 'Tail','Both',varargin);

if isempty(rperm) && strcmpi(Type,'splithalf')
    rperm=zeros(Reps,size(Data1,2));
    for i = 1:Reps
        rperm(i,:)=randperm(size(Data1,2));
    end
end

stats=struct;

if ~isempty(filterVal)
    Data1(Data1==filterVal)=nan;
    Data2(Data2==filterVal)=nan;
end
if strcmpi(Norm,'zscore')==1
    Data1=nan_zscore(Data1);
    Data2=nan_zscore(Data2);
end
if strcmpi(Norm,'fisherZ')==1
    Data1=atanh(Data1);
    Data2=atanh(Data2);
end

if strcmpi(Norm,'fisherZ&zscore')==1
    Data1=atanh(Data1);
    Data2=atanh(Data2);
    Data1=nan_zscore(Data1);
    Data2=nan_zscore(Data2);
end

if any(isnan(Data1(:))) || any(isnan(Data2(:)))
    if strcmpi(Type,'leave1out')
        rVals=nan(1,size(Data1,2));
        for i = 1:size(Data1,2)
            tData2=Data2;
            tData2(:,i)=[];
            if strcmp(corrType,'kendall')
                [rVals(1,i)]=rankCorr_Kendall_taua(Data1(:,i),nanmean(tData2,2));
            else
                [rVals(1,i)]=corr(Data1(:,i),nanmean(tData2,2),'type',corrType,'rows','pairwise');
            end    
        end         
    elseif strcmpi(Type,'splithalf')
        hlf=ceil(size(Data1,2)/2);
        rVals=nan(1,Reps);
        if strcmpi(corrType,'kendall')
            for i = 1:Reps
                [rVals(1,i)]=rankCorr_Kendall_taua(nanmean(Data1(:,rperm(i,1:hlf)),2),nanmean(Data2(:,rperm(i,hlf+1:end)),2));
            end
        else
            for i = 1:Reps
                [rVals(1,i)]=corr(nanmean(Data1(:,rperm(i,1:hlf)),2),nanmean(Data2(:,rperm(i,hlf+1:end)),2),'type',corrType,'rows','pairwise');
            end                
        end
    elseif strcmpi(Type,'leave1out_twoway') 
        rVals=nan(2,size(Data1,2));
        for i = 1:size(Data1,2)
            tData1=Data1;
            tData1(:,i)=[];
            tData2=Data2;
            tData2(:,i)=[];
            if strcmpi(corrType,'kendall')
                [rVals(1,i)]=rankCorr_Kendall_taua(Data1(:,i),nanmean(tData2,2));
                [rVals(2,i)]=rankCorr_Kendall_taua(Data2(:,i),nanmean(tData1,2));
            else
                [rVals(1,i)]=corr(Data1(:,i),nanmean(tData2,2),'type',corrType,'rows','pairwise');
                [rVals(2,i)]=corr(Data2(:,i),nanmean(tData1,2),'type',corrType,'rows','pairwise');
            end               
        end
        rVals=tanh(nanmean(atanh(rVals),1));
    end
else
    if strcmpi(Type,'leave1out')
        rVals=nan(1,size(Data1,2));
        for i = 1:size(Data1,2)
            tData2=Data2;
            tData2(:,i)=[];
            if strcmpi(corrType,'kendall')
                [rVals(1,i)]=rankCorr_Kendall_taua(Data1(:,i),mean(tData2,2));
            else
                [rVals(1,i)]=corr(Data1(:,i),mean(tData2,2),'type',corrType);                    
            end
        end  
    elseif strcmpi(Type,'splithalf')
        hlf=ceil(size(Data1,2)/2);
        rVals=nan(1,Reps);
        if strcmp(corrType,'kendall')
            for i = 1:Reps                
                [rVals(1,i)]=rankCorr_Kendall_taua(mean(Data1(:,rperm(i,1:hlf)),2),mean(Data2(:,rperm(i,hlf+1:end)),2));
            end     
        else
            for i = 1:Reps
                [rVals(1,i)]=corr(mean(Data1(:,rperm(i,1:hlf)),2),mean(Data2(:,rperm(i,hlf+1:end)),2),'type',corrType);
            end                 
        end  
    elseif strcmpi(Type,'leave1out_twoway') 
        rVals=nan(2,size(Data1,2));
        for i = 1:size(Data1,2)
            tData1=Data1;
            tData1(:,i)=[];
            tData2=Data2;
            tData2(:,i)=[];
            if strcmpi(corrType,'kendall')
                [rVals(1,i)]=rankCorr_Kendall_taua(Data1(:,i),mean(tData2,2));
                [rVals(2,i)]=rankCorr_Kendall_taua(Data2(:,i),mean(tData1,2));
            else
                [rVals(1,i)]=corr(Data1(:,i),mean(tData2,2),'type',corrType);  
                [rVals(2,i)]=corr(Data2(:,i),mean(tData1,2),'type',corrType); 
            end
        end 
        rVals=tanh(nanmean(atanh(rVals),1));
    end 
end    
     
rVals_clean=rVals;
rVals_clean(:,isnan(rVals_clean))=[];
rVals_clean(:,isinf(rVals_clean))=[];
if sum(contains(whichstats,'rmean'))==1
    stats.mean=mean(rVals_clean,2);
end
if sum(contains(whichstats,'fisherZmean'))==1
    stats.fisherZmean=mean(atanh(rVals_clean),2);
end

if sum(contains(whichstats,'tstat'))==1
    [~,stats.p,stats.CI,tstats]=ttest(atanh(rVals_clean'),0,'tail',Tail);
    stats.t=tstats.tstat;
    stats.sd=tstats.sd;
    stats.df=tstats.df;
end 

if sum(contains(whichstats,'prctiles'))==1
    prctiles=prctile(rVals_clean',[0,100,0.5,99.5,2.5,97.5,5,95,25,75,50]);
    stats.range=prctiles(1,[1,2]);
    stats.CI_99=prctiles(1,[3,4]);
    stats.CI_95=prctiles(1,[5,6]);
    stats.CI_90=prctiles(1,[7,8]);
    stats.CI_50=prctiles(1,[9,10]);
    stats.median=prctiles(1,11);
end 

if sum(contains(whichstats,'prctiles'))==1
    prctiles=prctile(atanh(rVals_clean)',[0,100,0.5,99.5,2.5,97.5,5,95,25,75,50]);
    stats.fzrange=prctiles(1,[1,2]);
    stats.fzCI_99=prctiles(1,[3,4]);
    stats.fzCI_95=prctiles(1,[5,6]);
    stats.fzCI_90=prctiles(1,[7,8]);
    stats.fzCI_50=prctiles(1,[9,10]);
    stats.fzmedian=prctiles(1,11);
end 