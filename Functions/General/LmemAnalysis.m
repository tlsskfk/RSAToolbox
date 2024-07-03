function [AllResultsTable] = LmemAnalysis(YTable,XTable,GroupingVar,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[figThresh] = VariableSetter('figThresh',0.05,varargin);
[SaveName] = VariableSetter('SaveName',[],varargin);
[SaveDir] = VariableSetter('SaveDir',[],varargin);
if isempty(SaveDir)
    SaveDir=uigetdir;
end
SaveDir=[SaveDir,'/'];
if ~exist(SaveDir)
    mkdir(SaveDir)
end 
if isempty(SaveName)
    SaveName=uiEnterName(['LMEM_Analysis_',genDateString],'Enter file name:');
end
YVarNames=YTable.Properties.VariableNames;
XVarNames=XTable.Properties.VariableNames;
if ~istable(GroupingVar)
    GroupingVar=cell2table(GroupingVar);
end
AllTables=cell(length(YVarNames)*length(XVarNames),1);
AllFormulas=cell(length(YVarNames)*length(XVarNames),1);
AllyNames=cell(length(YVarNames)*length(XVarNames),1);
AllxNames=cell(length(YVarNames)*length(XVarNames),1);
AllResults=cell(length(YVarNames)*length(XVarNames),1);

GroupNames=unique(GroupingVar.GroupingVar);
AllCorrs=nan(length(YVarNames)*length(XVarNames),length(GroupNames)*4);
GroupCorrNames=cell(1,length(GroupNames)*4);
count=1;
for i = 1:4:length(GroupNames)*4
    GroupCorrNames{1,i}=[GroupNames{count},'_r'];
    GroupCorrNames{1,i+1}=[GroupNames{count},'_p'];
    GroupCorrNames{1,i+2}=[GroupNames{count},'_s_r'];
    GroupCorrNames{1,i+3}=[GroupNames{count},'_s_p'];
    count=count+1;
end    
count=1;    
for YVarNum=1:length(YVarNames)
    yName=YVarNames{YVarNum};    
    for XVarNum=1:length(XVarNames)
        xName=XVarNames{XVarNum};
        AllyNames{count,1}=yName;
        AllxNames{count,1}=xName;
        AllFormulas{count,1}=[yName,' ~ ',xName,' + ','(1 + ',xName,'|GroupingVar)'];
        AllTables{count,1}=[YTable(:,yName),XTable(:,xName),GroupingVar];
        count2=1;
        for i = 1:4:length(GroupNames)*4
            YVar=YTable.(yName)(ismember(GroupingVar.GroupingVar,GroupNames{count2}),:);
            XVar=XTable.(xName)(ismember(GroupingVar.GroupingVar,GroupNames{count2}),:);
            [AllCorrs(count,i),AllCorrs(count,i+1)]=corr(YVar,XVar,'rows','pairwise');
            [AllCorrs(count,i+2),AllCorrs(count,i+3)]=corr(YVar,XVar,'rows','pairwise','type','spearman');
            count2=count2+1;
        end   
        count=count+1;
    end
end
corrTable=array2table(AllCorrs,'VariableNames',GroupCorrNames);
parfor i = 1:size(AllTables,1)
    [tempResults] = LME_PlusContinuousStats(AllTables{i,1},AllFormulas{i,1});
    AllResults{i,1}=tempResults(AllxNames{i,1},:);
end

for i = 1:size(AllResults,1)
    tempTable=AllResults{i,1};
    tempTable.Properties.RowNames={num2str(i)};
    tempTable=[cell2table({AllyNames{i,1},AllxNames{i,1}},'VariableNames',{'YName','XName'}),tempTable];
    if i == 1
        AllResultsTable = tempTable;
    else
        AllResultsTable=[AllResultsTable;tempTable];
    end
end
AllResultsTable=[AllResultsTable,corrTable];
figDir=[SaveDir,SaveName,'/'];
if ~exist(figDir)
    mkdir(figDir)
end 
for i = 1:height(AllResultsTable)
    if AllResultsTable.lme_pValue(i,1)<figThresh
        YName=AllResultsTable.YName{i,1};
        XName=AllResultsTable.XName{i,1};
        YVar=YTable.(YName);
        XVar=XTable.(XName);
        statsText=['Overall Stats',newline,' lme coef: ',num2str(AllResultsTable.lme_Estimate(i,1),4),newline,...
            ' lme tStat: ',num2str(AllResultsTable.lme_tStat(i,1),2),'; p =  ',num2str(AllResultsTable.lme_pValue(i,1),4),newline,...
            ' r = ',num2str(AllResultsTable.pearson_r(i,1),2),'; rho = ',num2str(AllResultsTable.spearman_r(i,1),2)];
        for j = 1:length(GroupNames)
            GroupName=GroupNames{j};
            statsText=[statsText,newline,GroupName,': r = ',num2str(AllResultsTable.([GroupName,'_r'])(i,1),2),'; p = ',num2str(AllResultsTable.([GroupName,'_p'])(i,1),2),newline,...
                ' rho = ',num2str(AllResultsTable.([GroupName,'_s_r'])(i,1),2),'; p = ',num2str(AllResultsTable.([GroupName,'_s_p'])(i,1),2)];
        end
        ScatterTitle='';
        [fs] = lmem_scatter(GroupingVar.GroupingVar,YVar,XVar,YName,XName,ScatterTitle,statsText);   
        figSaveName=[YName,'_',XName];
        set(gcf, 'Position', get(0, 'Screensize'));        
        export_fig([figDir,figSaveName,'.png'],'-Transparent','-png','-m2')
        close
    end
    
end
    
save([SaveDir,SaveName],'AllResultsTable','YTable','XTable','GroupingVar','AllFormulas');
end
    
