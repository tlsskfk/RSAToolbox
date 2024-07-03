function [fmriprep_table,AnalysisParameters] = Create_fMRIPrep_Table(ExperimentsDir,AnalysisParameters)
%Written by David Rothlein
%Input the directory of the fMRIPrep output. This function extracts
%information regarding particpants, tasks, runs filenames and directories etc...
%for the contained files. Output in a number of usable; however, only the 
%BIDsTable format is used for the subsequent further analyses.
%ExpDir Folder structure should be as follows: * = optional.
%ExpDir/[experiments]/fmriprep/[sub-ID]/[ses-0#]*/[anat,func]/files

%Logic of function: divide file properties into types,folders and properties,
%types consist of behavior, anat, or func, folders contain experiment, subID,
%sesion and file type. Properties consist of any infor contained in file
%name.
%All this info will be automatically compiled into a large table for each 
%file type and then relevant attribuites will be manually selected to construct
%the fMRIPrep_Table. This table will contain the local directory and relevant
%properties for each analysis unit (fMRI run).

%% select base directory containing Exp. folders(should be named Experiments)
if nargin==0
    ExperimentsDir=[];
    AnalysisParameters=struct;
elseif nargin==1
    AnalysisParameters=struct;
end
if isempty(ExperimentsDir)
    ExperimentsDir=uigetdir('/');
end


ssInfo=struct;
localDirs=[];
%% 
expDirInfo=dir(ExperimentsDir);
expDirInfo([1,2],:)=[];
for expNum = 1:size(expDirInfo,1)
    ssInfo.experiment={expDirInfo(expNum).name};
    tempDir=[ExperimentsDir,'/',ssInfo.experiment{1,1},'/fmriprep/'];
    if exist(tempDir)==0
        tempDir=[ExperimentsDir,'/',ssInfo.experiment{1,1},'/derivatives/fmriprep/'];
        fmriPrepFolder='/derivatives/fmriprep/';
        matlabFolder='/derivatives/matlab/';
    else
        fmriPrepFolder='/fmriprep/';
        matlabFolder='/matlab/';
    end
    if exist(tempDir)==0
        disp(['Unable to find fmriprep folder for experiment: ',expDirInfo(expNum).name, newline,...
        'Make sure format is = ExpDir/ExpName/fmriprep/... or ExpDir/ExpName/derivatives/fmriprep/...']);
        continue
    end
    subDirInfo=dir(tempDir);
    for ssNum = 1:size(subDirInfo,1)
        sesAnat=0; %Flag indicating whether anat data is session specific (0 = no) 
        if subDirInfo(ssNum).isdir==1 && contains(subDirInfo(ssNum).name,'sub')
            ssInfo.subName={subDirInfo(ssNum).name};
            ssInfo.sub={strrep(subDirInfo(ssNum).name,'sub-','')};
            ssInfo.UniqueID={[ssInfo.subName{1,1},'_',ssInfo.experiment{1,1}]};
            tempDir=[ExperimentsDir,'/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/'];
            sesDirInfo=dir(tempDir);
            for sesNum = 1:size(sesDirInfo,1)
                if sesDirInfo(sesNum).isdir==1 && contains(sesDirInfo(sesNum).name,'anat')
                    tempAnatDir=dir([ExperimentsDir,'/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/anat/']);
                    if size(tempAnatDir,1) > 7 %anat folder should have more than 7 files (tho cut off is arbitrary)
                        ssInfo.anatDir={['/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/anat/']};
                        sesAnat=1;
                    end
                    break
                end
            end
            for sesNum = 1:size(sesDirInfo,1)
                if sesDirInfo(sesNum).isdir==1 && contains(sesDirInfo(sesNum).name,'ses')
                    ssInfo.sesName={sesDirInfo(sesNum).name};
                    ssInfo.session=str2num(strrep(sesDirInfo(sesNum).name,'ses-',''));
                    if sesAnat == 0
                        ssInfo.anatDir={['/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/',ssInfo.sesName{1,1},'/anat/']};
                    end
                    ssInfo.funcDir={['/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/',ssInfo.sesName{1,1},'/func/']};
                    ssInfo.matDir={['/',ssInfo.experiment{1,1},matlabFolder,ssInfo.subName{1,1},'/',ssInfo.sesName{1,1}]};
                    localDirs=[localDirs;struct2table(ssInfo)];                    
                elseif sesDirInfo(sesNum).isdir==1 && contains(sesDirInfo(sesNum).name,'func')
                    ssInfo.sesName={'ses-01'};
                    ssInfo.session=1;
                    ssInfo.anatDir={['/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/anat/']};
                    ssInfo.funcDir={['/',ssInfo.experiment{1,1},fmriPrepFolder,ssInfo.subName{1,1},'/func/']};
                    ssInfo.matDir={['/',ssInfo.experiment{1,1},matlabFolder,ssInfo.subName{1,1},'/']};
                    localDirs=[localDirs;struct2table(ssInfo)];
                    break
                end
            end      
        end
    end
end

%% Select experiments, sessions, and participants to include in table
[FilterNames] = uiNameSelect(unique(localDirs.experiment),'Select experiments to include');
localDirs=localDirs(ismember(localDirs.experiment,FilterNames),:);
[FilterNames] = uiNameSelect(unique(localDirs.sesName),'Select sessions to include');
localDirs=localDirs(ismember(localDirs.sesName,FilterNames),:);
[FilterNames] = uiNameSelect(unique(localDirs.UniqueID),'Select participants to include');
localDirs=localDirs(ismember(localDirs.UniqueID,FilterNames),:);

%% Define which space analysis will be in (MNI, native, or other)
experiments=unique(localDirs.experiment);
if ~iscell(experiments)
    experiments={experiments};
end
experiments=experiments(:);
numExps=length(experiments);
count=1;
spaces=cell(1);
tasks=cell(1,3);
for expNum=1:numExps
    experiment=experiments{expNum,1};
    tempDirs=localDirs(ismember(localDirs.experiment,experiment),:);
    numSS=height(tempDirs);
    if numSS>5
        sampleSS=randperm(numSS);
        sampleSS=sampleSS(1,1:5); %select 5 random participants to identify space options
    else
        sampleSS=[1:numSS];
    end
    for ss=sampleSS
        funcFiles=dir([ExperimentsDir,tempDirs.funcDir{ss,1}]);
        for fileNum=1:size(funcFiles,1)
            if funcFiles(fileNum).isdir==0 && contains(funcFiles(fileNum).name,'space') && contains(funcFiles(fileNum).name,'task')
                try
                    tempInfo = fMRIPrep_Name2Info(funcFiles(fileNum).name);
                catch
                    continue
                end
                spaces{count,1}=strrep(tempInfo.space,'_boldref','');
                tasks{count,1}=['Exp- ',experiment,'; Task- ',tempInfo.task];
                tasks{count,2}=experiment;
                tasks{count,3}=tempInfo.task;
                count=count+1;
            end
        end
    end  
end
[~,taskInd]=unique(tasks(:,1));
tasks=tasks(taskInd,:);
singleSelect=1;
[UseSpace] = uiNameSelect(unique(spaces),'Select space to use',singleSelect);
singleSelect=0;
[UseTasks] = uiNameSelect(tasks(:,1),'Select task to use',singleSelect);
UseTasks=tasks(ismember(tasks(:,1),UseTasks),2:3);

%% Add default anat files (probseg masks) to localDirs table
default_anatFiles={'CSF_probseg','GM_probseg','WM_probseg'};
numSS=height(localDirs);
noHIT=1;
for SSnum = 1:numSS
    anatFiles=dir([ExperimentsDir,localDirs.anatDir{SSnum,1}]);
    for fileNum=1:size(anatFiles,1)
        if anatFiles(fileNum).isdir==0 && contains(anatFiles(fileNum).name,default_anatFiles)
            tempInfo = fMRIPrep_Name2Info(anatFiles(fileNum).name);
            if any(ismember('space', tempInfo.Properties.VariableNames))
                if strcmpi(tempInfo.space,UseSpace)
                    localDirs.(tempInfo.label){SSnum,1}=anatFiles(fileNum).name;
                    noHIT=0;
                end
            end
        end    
    end
end
if noHIT==1
    default_anatFiles={'label-CSF','label-GM','label-WM'};
    for SSnum = 1:numSS
        anatFiles=dir([ExperimentsDir,localDirs.anatDir{SSnum,1}]);
        for fileNum=1:size(anatFiles,1)
            if anatFiles(fileNum).isdir==0 && contains(anatFiles(fileNum).name,default_anatFiles)
                tempInfo = fMRIPrep_Name2Info(anatFiles(fileNum).name);
                if any(ismember('space', tempInfo.Properties.VariableNames))
                    if strcmpi(tempInfo.space,UseSpace)
                        localDirs.(tempInfo.label){SSnum,1}=anatFiles(fileNum).name;
                        noHIT=0;
                    end
                end
            end    
        end
    end
    default_anatFiles={'CSF','GM','WM'};
end

%% Create fmriprep_table by adding required func filenames for each exp, task, and run
default_funcFiles={'confounds_regressors','brain_mask','preproc_bold'};
default_funcFiles2={'confounds_timeseries','brain_mask','preproc_bold'};
fmriprep_table=table;
numTasks=size(UseTasks,1);
warning('off','MATLAB:table:RowsAddedExistingVars')
for taskNum = 1:numTasks
    experiment=UseTasks{taskNum,1};
    taskName=UseTasks{taskNum,2};
    tempDirs=localDirs(ismember(localDirs.experiment,experiment),:);    
    numSS=height(tempDirs);
    for SSnum = 1:numSS
        funcFiles=dir([ExperimentsDir,tempDirs.funcDir{SSnum,1}]);
        tempFileNames=table;
        fileFilter=[];
        for fileNum=1:size(funcFiles,1)
            if funcFiles(fileNum).isdir==0 && contains(funcFiles(fileNum).name,default_funcFiles) && contains(funcFiles(fileNum).name,['task-',taskName])
                tempInfo = fMRIPrep_Name2Info(funcFiles(fileNum).name);
                if ~strcmpi(tempInfo.task,taskName)
                    continue
                end
                %ensure pulled file names are in the selected space. if the
                %filename doesnt reference space, filename is pulled
                if any(ismember('space', tempInfo.Properties.VariableNames))
                    if ~strcmpi(tempInfo.space,UseSpace)
                        continue
                    end
                end
                %if filename doesn't indicate run number. assign run number
                %as 1. Otherwise pull run string and convert to number
                if ~any(ismember('run', tempInfo.Properties.VariableNames))
                    run=1;
                else
                    run=str2num(tempInfo.run);
                end
                tempFileNames.(tempInfo.desc){run,1}=funcFiles(fileNum).name;
                tempFileNames.run(run,1)=run;
                fileFilter(run,find(ismember(default_funcFiles,tempInfo.desc)))=1;
            end
        end
        fileFilter=sum(fileFilter,2)==length(default_funcFiles2);
        if sum(single(fileFilter(:)))==0
            funcFiles=dir([ExperimentsDir,tempDirs.funcDir{SSnum,1}]);
            tempFileNames=table;
            fileFilter=[];
            for fileNum=1:size(funcFiles,1)
                if funcFiles(fileNum).isdir==0 && contains(funcFiles(fileNum).name,default_funcFiles2) && contains(funcFiles(fileNum).name,['task-',taskName])
                    tempInfo = fMRIPrep_Name2Info(funcFiles(fileNum).name);
                    if ~strcmpi(tempInfo.task,taskName)
                        continue
                    end
                    %ensure pulled file names are in the selected space. if the
                    %filename doesnt reference space, filename is pulled
                    if any(ismember('space', tempInfo.Properties.VariableNames))
                        if ~strcmpi(tempInfo.space,UseSpace)
                            continue
                        end
                    end
                    %if filename doesn't indicate run number. assign run number
                    %as 1. Otherwise pull run string and convert to number
                    if ~any(ismember('run', tempInfo.Properties.VariableNames))
                        run=1;
                    else
                        run=str2num(tempInfo.run);
                    end
                    tempFileNames.(tempInfo.desc){run,1}=funcFiles(fileNum).name;
                    tempFileNames.run(run,1)=run;
                    fileFilter(run,find(ismember(default_funcFiles2,tempInfo.desc)))=1;
                end
            end
            fileFilter=sum(fileFilter,2)==length(default_funcFiles2);
        end
        tempFileNames=tempFileNames(fileFilter,:);
        numRuns=height(tempFileNames); 
        if numRuns==0
            continue
        end
        SS_Table=tempDirs(SSnum,[{'experiment'},{'sub'},{'session'},{'anatDir'},{'funcDir'},{'matDir'},default_anatFiles]);       
        SS_Table.task{1,1}=taskName;
        SS_Table.numRuns(1,1)=numRuns;
        SS_Table=[repmat(SS_Table,[numRuns,1]),tempFileNames];
%        try
            fmriprep_table=[fmriprep_table;SS_Table];
%         catch
%             fmriprep_table=SS_Table;
%         end
    end
end
VarNames=fmriprep_table.Properties.VariableNames;
if any(ismember(VarNames,{'GM','CSF','WM'}))
    VarNames=strrepCell(VarNames,'GM','GM_probseg');
    VarNames=strrepCell(VarNames,'WM','WM_probseg');
    VarNames=strrepCell(VarNames,'CSF','CSF_probseg');
end    
fmriprep_table.Properties.VariableNames=VarNames;
end