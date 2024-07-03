function [fmriprep_table] = gradCPT_BehaviorDataMigrate(fmriprep_table,ExperimentsDir,BehavDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Rename and move behavioral data
if nargin==2
    BehavDir=uigetdir('/');
    BehavDir=[BehavDir,'/'];
end    
  
a=ismember(fmriprep_table.experiment,'TRT_Reward') & ismember(fmriprep_table.task,'rcpt_final');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/TRT_Reward/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'Data_0_',tempBIDs.orig_IDs{i,1},'.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'


a=ismember(fmriprep_table.experiment,'TRACTS') & ismember(fmriprep_table.task,'gradCPT');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/TRACTS/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'Data_0_',tempBIDs.orig_IDs{i,1},'_v1B_fMRI_large.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'TRACTs') & ismember(fmriprep_table.task,'gradCPT');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/TRACTs/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'Data_0_',tempBIDs.orig_IDs{i,1},'_v1B_fMRI_large.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'OriginalGradCPT') & ismember(fmriprep_table.task,'gradCPT');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/OrigCC/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'Data_0_',tempBIDs.orig_IDs{i,1},'_',num2str(tempBIDs.run(i,1)),'.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'GradCPT_reward') & ismember(fmriprep_table.task,'gradCPTwithReward');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/Reward/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'Data_',tempBIDs.orig_IDs{i,1},'_',num2str(tempBIDs.run(i,1)),'.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'GradCPT_MindWandering') & ismember(fmriprep_table.task,'gradCPT');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/MW/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,tempBIDs.orig_IDs{i,1},'_Behavior_gradCPT_noTP.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'GradCPT_MindWandering') & ismember(fmriprep_table.task,'gradCPTMW');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/MW/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,tempBIDs.orig_IDs{i,1},'_Behavior_gradCPT_TP',num2str(tempBIDs.run(i,1)),'.mat'];
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end
    beh_raw=load(loadName);
    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'

a=ismember(fmriprep_table.experiment,'HalkoTMS') & ismember(fmriprep_table.task,'gradCPT');
tempBIDs=fmriprep_table(a,:);
loadbehDir=[BehavDir,'/HalkoTMS/'];
for i = 1:height(tempBIDs)
    loadName=[loadbehDir,'sub-',tempBIDs.sub{i,1},'_ses-',num2str(tempBIDs.session(i,1)),'_run-',num2str(tempBIDs.run(i,1)),'.mat'];    
    try
        beh_raw=load(loadName);
    catch
        disp(['No file: ',loadName])
        continue
    end
    behDir=strrep(tempBIDs.funcDir{i,1},'/fmriprep/','/matlab/');
    behDir=strrep(behDir,'/func/',['/beh/raw/']);
    behDir=[ExperimentsDir,behDir];
    saveName=['sub-',tempBIDs.sub{i,1},'_task-',tempBIDs.task{i,1},'_run-',num2str(tempBIDs.run(i,1)),'_desc-beh_raw.mat'];
    if ~exist(behDir,'file')
        mkdir(behDir);
    elseif exist([behDir,saveName],'file')
        continue
    end

    save([behDir,saveName],'beh_raw');
end
clearvars -except 'fmriprep_table' 'BehavDir' 'ExperimentsDir'
end