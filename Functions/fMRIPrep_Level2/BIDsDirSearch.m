function [filePaths,AnalysisType,AnalysisName,ExperimentsDir,fmriprep_table,fmriprep_table_name,SubjectOrRun,AllAnalysisNames,AllAnalysisTypes] = BIDsDirSearch(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Pulls information about analyse names and generates file names and
%locations directly from the folder structure. Can also be used to delete
%or rename files.
%[filePaths]=BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',AnalysisType,'AnalysisName',AnalysisName,'SubjectOrRun',SubjectOrRun)


if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end

[AnalysisType] = VariableSetter('AnalysisType',[],varargin);
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
[GetAnalysisNames] = VariableSetter('GetAnalysisNames',0,varargin);
[GetAnalysisTypes] = VariableSetter('GetAnalysisTypes',0,varargin);
[AnalysisNameFilter] = VariableSetter('AnalysisNameFilter',[],varargin);
[DeleteFiles] = VariableSetter('DeleteFiles',0,varargin);
[DeleteAnalysisName] = VariableSetter('DeleteAnalysisName',0,varargin);
%[RenameFiles] = VariableSetter('RenameFiles',0,varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[TitleTextType] = VariableSetter('TitleTextType','Select type of analysis to search',varargin);
[TitleTextName] = VariableSetter('TitleTextName','Select name of analysis to compile',varargin);
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
filePaths=table;
warning('off','MATLAB:table:RowsAddedExistingVars');
tableVars=fmriprep_table.Properties.VariableNames;
SubOrRunOpts={'Subject','Run'};
if any(ismember(tableVars,'numRuns_bySes'))
    SubOrRunOpts=[SubOrRunOpts,{'Session'}];
end
if any(ismember(tableVars,'numRuns_byGroup'))
    SubOrRunOpts=[SubOrRunOpts,{'Group'}];
end
if ~any(ismember(tableVars,'numRuns_bySub'))
    fmriprep_table.numRuns_bySub=fmriprep_table.numRuns;
end
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect(SubOrRunOpts,'Perform analysis by subject or by run:',SingleSelect);
end
bySes=0;
bySub=0;
byGroup=0;
byRun=0;
groupName=[];
ssAppend=[];
if strcmpi(SubjectOrRun,'Subject')
    bySub=1;
    useIndicies=find(fmriprep_table.numRuns_bySub)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    ssAppend='_bySub';
elseif strcmpi(SubjectOrRun,'Group')
    byGroup=1;
    useIndicies=find(fmriprep_table.numRuns_byGroup)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_byGroup; 
    groupName=uiEnterName('','Enter group name.');
    ssAppend=['_by',groupName];
elseif strcmpi(SubjectOrRun,'Session')
    bySes=1;
    useIndicies=find(fmriprep_table.numRuns_bySes)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySes; 
    ssAppend='_bySes';
else    
    byRun=1;
    useIndicies=[1:TotalRuns];
    ssApend='_byRun';
end
AllAnalysisTypes=[];

if isempty(AnalysisType)
    SingleSelect=1;
    for i = useIndicies
        dirName=[ExperimentsDir,fmriprep_table.matDir{i,1}];
        [folderNames]=getFolderAndFileNames(dirName);
        AllAnalysisTypes=unique([AllAnalysisTypes;folderNames]);
    end
    inNames=AllAnalysisTypes;
    if GetAnalysisTypes == 0
        [AnalysisType] = uiNameSelect(inNames,TitleTextType,SingleSelect);
    else
        filePaths=AllAnalysisTypes;
    end
end

if iscell(AnalysisType)
    AnalysisType=AnalysisType{1,1};
end
AllAnalysisNames=[];
if GetAnalysisTypes == 0
if isempty(AnalysisName)
    SingleSelect=1;
    AllAnalysisNames=[];
    for i = useIndicies
        dirName=[ExperimentsDir,fmriprep_table.matDir{i,1},'/',AnalysisType,'/'];
        [folderNames]=getFolderAndFileNames(dirName);
        AllAnalysisNames=unique([AllAnalysisNames;folderNames]);
    end
    if ~isempty(AnalysisNameFilter)
        inNames=AllAnalysisNames(contains(AllAnalysisNames,AnalysisNameFilter),:);
    else
        inNames=AllAnalysisNames;
    end
    if GetAnalysisNames == 0
        [AnalysisName] = uiNameSelect(inNames,TitleTextName,SingleSelect);
    else
        filePaths=AllAnalysisNames;
    end
end

if iscell(AnalysisName)
    AnalysisName=AnalysisName{1,1};
end

if GetAnalysisNames == 0
for i = useIndicies
    dirName=[ExperimentsDir,fmriprep_table.matDir{i,1},'/',AnalysisType,'/',AnalysisName,'/'];
    if DeleteAnalysisName == 1
        try
            rmdir(dirName,'s');
            continue
        catch
            continue
        end
    end
    [parcelNames]=getFolderAndFileNames(dirName);
    
    if isempty(parcelNames)
        parcelNames={''};
        tableLabels={AnalysisName};
    elseif ~iscell(parcelNames)
        parcelNames={parcelNames}; 
        tableLabels=parcelNames;
    else
        parcelNames=parcelNames(:); 
        tableLabels=parcelNames;
    end

    for parcelNum=1:length(parcelNames)
        parcelName=parcelNames{parcelNum,1};
        tableLabel=tableLabels{parcelNum,1};
        dirName=[ExperimentsDir,fmriprep_table.matDir{i,1},'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
        dirName=strrep(dirName,'\','/');
        dirName=strrep(dirName,'//','/');
        dirName=strrep(dirName,'//','/');
        [~,fileNames]=getFolderAndFileNames(dirName);        
        if isempty(fileNames)
            if DeleteFiles==1
                delete(dirName);
            end
            continue
        end
        numRuns=fmriprep_table.numRuns(i,1);
        runNum=fmriprep_table.run(i,1);
        if numRuns==1
            for j = 1:length(fileNames)
                [fMRIPrepInfo] = fMRIPrep_Name2Info(fileNames{j,1});
                if any(ismember(fMRIPrepInfo.Properties.VariableNames,'run'))
                    fMRIPrepInfo(:,{'run'})=[];
                end
                MatchVars=fmriprep_table(i,fmriprep_table.Properties.VariableNames(ismember(fmriprep_table.Properties.VariableNames,fMRIPrepInfo.Properties.VariableNames)));
                MatchNames=fMRIPrepInfo(1,fMRIPrepInfo.Properties.VariableNames(ismember(fMRIPrepInfo.Properties.VariableNames,MatchVars.Properties.VariableNames)));
                if any(~ismember(table2cell(MatchVars),table2cell(MatchNames)))
                    continue
                else
                    filePaths.(tableLabel){i,1}=[dirName,fileNames{j,1}];
                    break
                end
            end                        
        elseif bySub==1
            for j = 1:length(fileNames)        
                [fMRIPrepInfo] = fMRIPrep_Name2Info(fileNames{j,1});
                if ~any(ismember(fMRIPrepInfo.Properties.VariableNames,'run'))
                    MatchVars=fmriprep_table(i,fmriprep_table.Properties.VariableNames(ismember(fmriprep_table.Properties.VariableNames,fMRIPrepInfo.Properties.VariableNames)));
                    MatchNames=fMRIPrepInfo(1,fMRIPrepInfo.Properties.VariableNames(ismember(fMRIPrepInfo.Properties.VariableNames,MatchVars.Properties.VariableNames)));
                    if any(~ismember(table2cell(MatchVars),table2cell(MatchNames)))
                        continue
                    else
                        if length(tableLabel)>namelengthmax
                            tableLabel=tableLabel(1,1:namelengthmax);
                        end
                        filePaths.(tableLabel){i,1}=[dirName,fileNames{j,1}];
                        break
                    end
                end
            end
        else
            for j = 1:length(fileNames)        
                [fMRIPrepInfo] = fMRIPrep_Name2Info(fileNames{j,1});
                if any(ismember(fMRIPrepInfo.Properties.VariableNames,'run')) && any(ismember(fMRIPrepInfo.Properties.VariableNames,'task'))
                    fileRunNum=str2num(fMRIPrepInfo.run);
                    if runNum==fileRunNum && strcmpi(fMRIPrepInfo.task,fmriprep_table.task{i,1})
                        filePaths.(tableLabel){i,1}=[dirName,fileNames{j,1}];
                        break
                    end
                end
            end        
        end
        if DeleteFiles==1
        	delete(filePaths.(tableLabel){i,1});
        end  
    end
end
end        
end
end

