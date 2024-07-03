function [fmriprep_table] = ComputeBehavioralVars(fmriprep_table,ExperimentsDir,varargin)
%Template function for data processing from the BIDsTable

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);

%For other variable inputs XX use VariableSetter function where 'Variable' is the
%text string indicating the variable name and DefaultVal is the value
%assigned to Variable if nothing is specified.
[RunSubSample] = VariableSetter('RunSubSample',0,varargin);


%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='beh_vars';

%Allows you to set name for this particular analysis
if RunSubSample==1
    if isempty(AnalysisName)
        AnalysisName=uiEnterName(['SubSample_'],['Enter name for ',AnalysisType,newline,'analysis below:']);
    end
else
   if isempty(AnalysisName)
       AnalysisName='Overall';
   end
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
    if length(unique(fmriprep_table.session))>1
        useIndicies=find(fmriprep_table.numRuns_bySes)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;        
    else
        useIndicies=find(fmriprep_table.numRuns_bySub)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    end
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end
%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%% Search a mat file for names to select from
%Alternatively use this function:
%[outNames] = BIDsDirSearch(fmriprep_table,ExperimentsDir,varargin);

%Select a subset of items from a larger list
%[EventNames] = uiNameSelect([UniqueEventNames],'Select events or blocks to include:');

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
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
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=['sub-',fmriprep_table.sub{dataInd,1},'_task-',fmriprep_table.task{dataInd,1},'_desc-beh_vars.mat'];
    else
        SaveName=['sub-',fmriprep_table.sub{dataInd,1},'_task-',fmriprep_table.task{dataInd,1},'_run-',num2str(fmriprep_table.run(dataInd,1)),'_desc-beh_vars.mat'];
    end
    %% If analysis does not involves parcellations:
    SavePrefix=[ExperimentsDir,SaveDir,'/'];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
    end    
    SaveNames{1,1}=[SavePrefix,SaveName];
    if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
        disp(['Skipping-- file exists: ',SaveName]);
        continue
    end    
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    InputData=cell(1);
    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        LoadVars={'beh_raw'};
        loadInd=dataInd+run-1;
        %skip previous errors
%         if ~isempty(fmriprep_table.Error{loadInd,1})
%             continue
%         end    
        %% set load paths and variable names
        LoadVars={'beh_raw'};
        LoadDir=strrep(fmriprep_table.funcDir{loadInd,1},'/fmriprep/','/matlab/');
        LoadDir=strrep(LoadDir,'/func/','/beh/raw/');
        LoadName=['sub-',fmriprep_table.sub{loadInd,1},'_task-',fmriprep_table.task{loadInd,1},'_run-',num2str(fmriprep_table.run(loadInd,1)),'_desc-beh_raw.mat'];
        LoadPath=[ExperimentsDir,LoadDir,LoadName];
        try
            TempLoadData = load(LoadPath,LoadVars{1,1});
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,LoadPath]);
            continue
        end       
        InputData{count,1}=LoadPath; 
        count=count+1;
    end  
    
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ']);
        continue
    end  
    
    %% Run analysis here!!
    try
        [ inFilesBySS ] = gradCPT_LoadBehavior(InputData);
        [ Behavior ] = gradCPT_SplitHalfBehavior( [],[],inFilesBySS{1,1},1,[],0);
        beh_vars=array2table(Behavior.FullOutput,'VariableNames',Behavior.FullOutputLabels);
        save(SaveNames{1,1},'beh_vars');    
    catch
        disp('Error!')
    end
    toc
end
end 

function [ inFilesBySS,inFilesByRun ] = gradCPT_LoadBehavior(loadNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


numSS=1;
numRuns=size(loadNames,1);
count=1;
for n=1:numSS 
    inFilesBySS{n,1}.VTC=[];
    inFilesBySS{n,1}.ttt=[];
    inFilesBySS{n,1}.response=[];
    inFilesBySS{n,1}.data=[];
    inFilesBySS{n,1}.Zone=[];
    inFilesBySS{n,1}.ZonePrime=[];
    inFilesBySS{n,1}.VTCderiv=[];
    inFilesBySS{n,1}.bordertracker=[];
    inFilesBySS{n,1}.numTrials=0;
    for r = 1:numRuns
        try
            load(loadNames{r,1}); %beh_raw
            inFilesByRun{1,count}.ttt=beh_raw.ttt;
            inFilesByRun{1,count}.response=beh_raw.response;
            inFilesByRun{1,count}.data=beh_raw.data; 
            if isfield(beh_raw,'bordertracker')
                inFilesByRun{1,count}.bordertracker=beh_raw.bordertracker;
            end
        catch
            continue
        end
        inFilesByRun{1,count}.VTC=CPT_analyze_zone_func2(inFilesByRun{1,count}.response,inFilesByRun{1,count}.data);
        inFilesByRun{1,count}.ttt(end,:)=[];
        inFilesByRun{1,count}.response(end,:)=[];
        inFilesByRun{1,count}.data(end,:)=[];
        inFilesByRun{1,count}.numTrials=length(inFilesByRun{1,count}.ttt);
        [inFilesByRun{1,count}.Zone,~,inFilesByRun{1,count}.ZonePrime,~,inFilesByRun{1,count}.VTCderiv] = ZoneByItem( inFilesByRun{1,count}.VTC,inFilesByRun{1,count}.data(:,3),inFilesByRun{1,count}.data(:,4));
        if isfield(inFilesByRun{1,count},'bordertracker')
            inFilesByRun{1,count}.bordertracker(end,:)=[];
            inFilesByRun{1,count}.ComputeReward=1;
        else
            inFilesByRun{1,count}.ComputeReward=0;
            inFilesByRun{1,count}.bordertracker=ones(inFilesByRun{1,count}.numTrials,3)*250;
        end
        if size(inFilesByRun{1,count}.data,2)<18
            inFilesByRun{1,count}.data=[inFilesByRun{1,count}.data,zeros(size(inFilesByRun{1,count}.data,1),18-size(inFilesByRun{1,count}.data,2))];
        end   
        inFilesBySS{n,1}.VTC=[inFilesBySS{n,1}.VTC;inFilesByRun{1,count}.VTC];
        inFilesBySS{n,1}.Zone=[inFilesBySS{n,1}.Zone;inFilesByRun{1,count}.Zone];
        inFilesBySS{n,1}.ZonePrime=[inFilesBySS{n,1}.ZonePrime;inFilesByRun{1,count}.ZonePrime];
        inFilesBySS{n,1}.VTCderiv=[inFilesBySS{n,1}.VTCderiv;inFilesByRun{1,count}.VTCderiv];        
        inFilesBySS{n,1}.ttt=[inFilesBySS{n,1}.ttt;inFilesByRun{1,count}.ttt];
        inFilesBySS{n,1}.response=[inFilesBySS{n,1}.response;inFilesByRun{1,count}.response];
        try
            inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;inFilesByRun{1,count}.data];
        catch
            if size(inFilesBySS{n,1}.data,2)<size(inFilesByRun{1,count}.data,2)
                inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;inFilesByRun{1,count}.data(:,1:size(inFilesBySS{n,1}.data,2))];
            else
                temp=zeros(size(inFilesByRun{1,count}.data,1),size(inFilesBySS{n,1}.data,2));
                temp(:,1:size(inFilesByRun{1,count}.data,2))=inFilesByRun{1,count}.data;
                inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;temp];
            end    
        end    
        inFilesBySS{n,1}.bordertracker=[inFilesBySS{n,1}.bordertracker;inFilesByRun{1,count}.bordertracker];  
        inFilesBySS{n,1}.ComputeReward=inFilesByRun{1,count}.ComputeReward;
        inFilesBySS{n,1}.numTrials=inFilesBySS{n,1}.numTrials+inFilesByRun{1,count}.numTrials;
        count=count+1;
    end
end
end