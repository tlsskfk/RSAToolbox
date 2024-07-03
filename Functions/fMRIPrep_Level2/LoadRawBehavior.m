function [fmriprep_table] = LoadRawBehavior(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Load behavior saved in multiple formats

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Run',varargin);
%Subject or run level analysis. Will prompt request.
[DataFormat] = VariableSetter('DataFormat',[],varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',['raw'],varargin);
%Remove first N rows of datafile 
[RemoveRows] = VariableSetter('RemoveRows',[],varargin);
%Variable names to save (cell of strings)
[SaveVariableNames] = VariableSetter('SaveVariableNames',[],varargin);
%Variable names to save (cell of strings)
[UseDefault] = VariableSetter('UseDefault',[],varargin);
[BaseFolder] = VariableSetter('BaseFolder',[],varargin);
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='beh';

%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName('raw',['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=dataInd;
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end

if isempty(DataFormat)
    SingleSelect=1; %Allows only a single value to be selected.
    [DataFormat] = uiNameSelect({'TrialType_And_Dur','EventTimecourse','AFNI','PsychoPy','CSV','TSV'},'Select Raw Data Format:',SingleSelect);
end
loadSuffix=uiEnterName('e.g. desc-behavraw or events','Enter rawfile suffix');
if strcmpi(DataFormat,'PsychoPy')|| strcmpi(DataFormat,'CSV')
     loadSuffix=[loadSuffix,'.csv'];
elseif strcmpi(DataFormat,'tsv')
     loadSuffix=[loadSuffix,'.tsv'];
end

if isempty(UseDefault)
    [UseDefault] = uiNameSelect({'Yes','No'},'Use default file name? (No to select)',1);
    if strcmpi(UseDefault,'Yes')
        UseDefault=1;
    else
        UseDefault=0;
    end
end
if UseDefault==1
    if isempty(BaseFolder)
        [BaseFolder] = uiNameSelect({'fmriprep','sourcedata','other'},'Select default behavior folder.',1);
        if strcmpi(BaseFolder,'other')
            BaseFolder='fmriprep';
            UseDefault=0;
        end
    end 
end
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
datadirSwitch=0;
if any(ismember(fmriprep_table.Properties.VariableNames,'confounds_regressors'))
    confoundName='confounds_regressors';
else
    confoundName='confounds_timeseries';
end
    
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    
    %% Set save directory and save name
    SaveDir=strrep(fmriprep_table.funcDir{dataInd,1},'/func/',['/',AnalysisType,'/',AnalysisName,'/']);
    SaveDir=strrep(SaveDir,'/fmriprep/','/matlab/');
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-beh_raw'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis does not involves parcellations:
    SavePrefix=[ExperimentsDir,SaveDir];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
    end    
    SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
        disp(['Skipping-- file exists: ',SaveName]);
        continue
    end    
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.NumRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
%         if ~isempty(fmriprep_table.Error{loadInd,1})
%             continue
%         end    
        %% set load paths and variable names        
        LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.funcDir{loadInd,1},strrep(fmriprep_table.(confoundName){loadInd,1},['desc-',confoundName,'.tsv'],loadSuffix)];
        LoadPaths{1,1}=strrep(LoadPaths{1,1},'fmriprep',BaseFolder);
        count=count+1;
    end  
    
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ']);
        continue
    end  
    
    %% Run analysis here!!
    if strcmpi(DataFormat,'PsychoPy')|| strcmpi(DataFormat,'CSV') || strcmpi(DataFormat,'TSV')
        
        try
            if strcmpi(DataFormat,'PsychoPy')|| strcmpi(DataFormat,'CSV')
                try
                    beh_raw=readtable(LoadPaths{1,1});
                catch
                    LoadPaths{1,1}=strrep(LoadPaths{1,1},'_run-','_run-0');
                    beh_raw=readtable(LoadPaths{1,1});
                end
            elseif strcmpi(DataFormat,'TSV')
                try
                    beh_raw=struct2table(tdfread(LoadPaths{1,1},'\t'));
                catch
                    LoadPaths{1,1}=strrep(LoadPaths{1,1},'_run-','_run-0');
                    beh_raw=struct2table(tdfread(LoadPaths{1,1},'\t'));
                end
            end
        catch
            if UseDefault==0
                if datadirSwitch==0
                    BehavDataDir=uigetdir(ExperimentsDir);
                    datadirSwitch=1;
                end
                if ismember(fmriprep_table.Properties.VariableNames,'ses')
                    disp([fmriprep_table.experiment{dataInd,1},' sub-',fmriprep_table.sub{dataInd,1},' ses-',fmriprep_table.ses{dataInd,1},' task-',fmriprep_table.task{dataInd,1},' run-',fmriprep_table.run{dataInd,1}]);
                else
                    disp([fmriprep_table.experiment{dataInd,1},' sub-',fmriprep_table.sub{dataInd,1},' task-',fmriprep_table.task{dataInd,1},' run-',fmriprep_table.run{dataInd,1}]);
                end
                [tempFileName,tempFilePath]=uigetfile([BehavDataDir,'*.*'],'Select behavioral data file');
                LoadPaths{1,1}=[tempFilePath,tempFileName];
                try
                    if strcmpi(DataFormat,'PsychoPy')|| strcmpi(DataFormat,'CSV')
                        beh_raw=readtable(LoadPaths{1,1});
                    elseif strcmpi(DataFormat,'TSV')
                        beh_raw=struct2table(tdfread(LoadPaths{1,1},'\t'));
                    end
                catch
                    disp(['Skipping subject-- data will not load!']);
                    continue
                end
            else
                disp(['Skipping subject-- data will not load!']);
                continue        
            end               
        end    
        if isempty(RemoveRows)
            RemoveRows=str2num(uiEnterName('0','Select number of rows to remove.')); 
        end
        if RemoveRows > 0
            beh_raw([1:RemoveRows],:)=[]; 
        end        
        if isempty(SaveVariableNames)
            SaveVariableNames=uiNameSelect(beh_raw.Properties.VariableNames);
        end
        beh_raw=beh_raw(:,SaveVariableNames);
        SaveVariableNames=SaveVariableNames(:);
        for j=1:size(SaveVariableNames,1)
            if ischar(beh_raw.(SaveVariableNames{j,1}))
                beh_raw.(SaveVariableNames{j,1})=Char2Cell(beh_raw.(SaveVariableNames{j,1}));
            end
        end
        save(SaveNames{1,1},'beh_raw');
    end     
    toc
end
end

function NewCell=Char2Cell(OldChar)
    numItems=size(OldChar,1);
    NewCell=cell(numItems,1);
    for i=1:numItems
        tempStr=OldChar(i,:);
        tempStr(strfind(tempStr,' '))=[];
        NewCell{i,1}=tempStr;
    end
end
 
