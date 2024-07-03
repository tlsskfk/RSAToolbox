function [behRSMs_NN,behRSMs_Mean,behRSMs_Min,behRSMs_Max,behRSMs_absDiffMean,behRSMs_absDiffMeanInv] = ssRSMs_BehavioralVars(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%Compile Behavioral Data

[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Subject',varargin);
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
[behVars_names] = VariableSetter('behVars_names',[],varargin);
[subGroups] = VariableSetter('subGroups',[],varargin);
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
GroupDir=strrep(ExperimentsDir,'Experiments/',['GroupAnalysis_ssRSMs/',fmriprep_table_name,'/Behavior/']);
GroupDir=strrep(GroupDir,'//','/');
if ~exist(GroupDir,'file')
    mkdir(GroupDir);     
end

try
    [taskVars,~,taskVarInd] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','beh_vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','beh_vars',...
        'AnalysisName','Overall',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    IncludeTaskVars=1;
catch
    IncludeTaskVars=0;
end

try
    [ssVars,~,ssVarInd] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','SS_Vars',...
        'LoadVarFormat','Table',...
        'DataType','Other',...
        'AnalysisType','SSVars',...
        'AnalysisName','SSVars_All',...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Table',...
        'ParcelName','None');
    IncludeSSVars=1;
catch
    IncludeSSVars=0;
end

if IncludeSSVars + IncludeTaskVars == 2
    if sum(ssVarInd) > sum(taskVarInd)
        ssVars=ssVars(taskVarInd==1,:);
        used_fmriprep_table=fmriprep_table(taskVarInd==1,:);
        selectInd=taskVarInd;
    elseif sum(ssVarInd) < sum(taskVarInd)
        taskVars=taskVars(ssVarInd==1,:);
        used_fmriprep_table=fmriprep_table(ssVarInd==1,:);
        selectInd=ssVarInd;
    end   
    behVars=[taskVars,ssVars];
elseif IncludeSSVars == 1
    selectInd=ssVarInd;
    behVars=ssVars;
elseif IncludeTaskVars == 1
    selectInd=taskVarInd;    
    behVars=taskVars;
else
    disp('No behavioral data!!!!')
    return
end
if isempty(behVars_names)
    behVars_names=behVars.Properties.VariableNames;
    behVars_names=uiNameSelect(behVars_names,'Select behavioral varaibles to compare.',0);
    removeInd=zeros(length(behVars_names),1);
    N=height(behVars);
    for i = 1:length(behVars_names)
        tempVar=behVars.(behVars_names{i,1});
        if ~isnumeric(tempVar)
            removeInd(i,1)=1;
            continue
        end
        if all(isnan(tempVar))
            removeInd(i,1)=1;
            continue
        end    
        if all(isinf(tempVar))
            removeInd(i,1)=1;
            continue
        end  
        [~,F] = mode(tempVar);
        if F == N
            removeInd(i,1)=1;
            continue
        end
    end
    behVars_names(removeInd==1,:)=[];
end
behVars=behVars(:,behVars_names);
if isempty(subGroups)
    subGroupNum=uiEnterName('0','Enter number of subgroups to compile');
    subGroupNum=str2num(subGroupNum);
    if subGroupNum ~= 0
        subGroups=cell(subGroupNum,2);
        for i = 1:subGroupNum
            subGroups{i,1}=uiNameSelect(behVars_names,'Select subgroup set.',0);
            subGroups{i,2}=uiEnterName('','Enter sub group name:');
        end
    else
        subGroups=[];
    end
else
   subGroupNum=size(subGroups,1);
end
[behRSMs_NN,behRSMs_Mean,behRSMs_Min,behRSMs_Max,behRSMs_absDiffMean,behRSMs_absDiffMeanInv]=Compute1D_RSMs(table2array(behVars),used_fmriprep_table.sub);
behRSMs_NN=array2table(behRSMs_NN,'VariableNames',behVars_names);
behRSMs_Mean=array2table(behRSMs_Mean,'VariableNames',behVars_names);
behRSMs_Min=array2table(behRSMs_Min,'VariableNames',behVars_names);
behRSMs_Max=array2table(behRSMs_Max,'VariableNames',behVars_names);
behRSMs_absDiffMean=array2table(behRSMs_absDiffMean,'VariableNames',behVars_names);
behRSMs_absDiffMeanInv=array2table(behRSMs_absDiffMeanInv,'VariableNames',behVars_names);
save([GroupDir,'AllVars'],'behRSMs_NN','behRSMs_Mean','behRSMs_Min','behRSMs_Max','behRSMs_absDiffMean','behRSMs_absDiffMeanInv');
save([GroupDir,'used_fmriprep_table'],'used_fmriprep_table','selectInd');
VarGroupDir = [GroupDir,'byVar/'];
if ~exist(VarGroupDir,'file')
    mkdir(VarGroupDir);     
end
RSMTypeNames={'NN','Mean','Min','Max','absDiffMean','behRSMs_absDiffMeanInv'};
for i = 1:length(behVars_names)
    varName=behVars_names{i,1};
    tempArray=[behRSMs_NN.(varName),behRSMs_Mean.(varName),behRSMs_Min.(varName),behRSMs_absDiffMean.(varName),behRSMs_Max.(varName),behRSMs_absDiffMeanInv.(varName)];
    behRSMs=array2table(tempArray,'VariableNames',RSMTypeNames);
    save([VarGroupDir,varName],'behRSMs','selectInd');
end

if ~isempty(subGroups)
    SubGroupDir = [GroupDir,'VarSets/'];
    if ~exist(SubGroupDir,'file')
        mkdir(SubGroupDir);     
    end  
    for i = 1:subGroupNum
        UsedVarNames=subGroups{i,1};
        UsedBehVars=table2array(behVars(:,UsedVarNames));
        VarSetRSMs=array2table(corr(UsedBehVars,'rows','pairwise'),'RowNames',UsedVarNames,'VariableNames',UsedVarNames);
        ssRSM_BehVarSet=mat2uppertriuvectormat(corr(scaleVals(UsedBehVars,1,0,1)','rows','pairwise'));        
        save([SubGroupDir,subGroups{i,2}],'VarSetRSMs','ssRSM_BehVarSet','selectInd','UsedVarNames','UsedBehVars');
    end    
end    
end
