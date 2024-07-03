function [fmriprep_table] = ComputeSubSampleBehavioralVars(fmriprep_table,ExperimentsDir,varargin)
%Template function for data processing from the BIDsTable

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
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
AnalysisType='subsample_beh_vars';

%Allows you to set name for this particular analysis
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
%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%% Search a mat file for names to select from
%Alternatively use this function:
%[outNames] = BIDsDirSearch(fmriprep_table,ExperimentsDir,varargin);
[filePaths_subsamples,~,AnalysisName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','subsample','TitleTextName',['Select subsamples for beh vars:']);
filePaths_subsamples=filePaths_subsamples.(AnalysisName);
[filePaths_behData] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','raw');
filePaths_behData=filePaths_behData.raw;
[save_filePaths] = strrepCell(filePaths_subsamples,'/subsample/',['/',AnalysisType,'/']);
[save_filePaths] = strrepCell(save_filePaths,'desc-subsample.',['desc-subsample_behvars.']);
if bySS==1
    [save_filePaths] = strrepCell(save_filePaths,'_run-01_',['_']);
end
%Select a subset of items from a larger list
%[EventNames] = uiNameSelect([UniqueEventNames],'Select events or blocks to include:');
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    SaveName=save_filePaths{dataInd,1};
    if bySS==1
        startRun=fmriprep_table.run(dataInd,1);
        if startRun < 10
            startRun=['0',num2str(startRun)];
        else
            startRun=num2str(startRun);
        end
        SaveName=strrep(SaveName,['_run-',startRun],'');
    end
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    [SaveDir] = filepath2DirAndFileName({SaveName});
    %% Set save directory and save name
    if ~exist(SaveDir{1,1},'file')
        mkdir(SaveDir{1,1});
    end    
    
    if exist(SaveName,'file')~=0 && Overwrite==0
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
    subsampleData=[];
    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        LoadVars{1,1}='SubAssignsByTrial';
        %skip previous errors
        if ~isempty(fmriprep_table.Error{loadInd,1})
            continue
        end    
        %% set load paths and variable names
        try
            TempLoadData = load(filePaths_subsamples{loadInd,1},LoadVars{1,1});
            TempLoadData=TempLoadData.(LoadVars{1,1});
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,filePaths_subsamples{loadInd,1}]);
            continue
        end     
        try
            TempBehData = load(filePaths_behData{loadInd,1},'beh_raw');
            TempBehData = TempBehData.beh_raw;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,filePaths_behData{loadInd,1}]);
            continue
        end  
        TempLoadData(TempLoadData==0)=nan;
        TempLoadData=fillmissing(TempLoadData,'nearest',1);
        subsampleData=[subsampleData;quikCellFormat(TempLoadData)];
        
        InputData{count,1}=filePaths_behData{loadInd,1}; 
        count=count+1;
    end  
    numSubs=size(subsampleData,2);
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ']);
        continue
    end  
    
    %% Run analysis here!!
    beh_vars = cell(numSubs,1);
    labels=cell(numSubs,1);
    tic
    parfor_progress(numSubs);
    parfor sub = 1:numSubs       
        [ inFilesBySS ] = gradCPT_LoadBehavior(InputData,subsampleData(:,sub));
        [ Behavior ] = gradCPT_SplitHalfBehavior( [],[],inFilesBySS{1,1},1,[],0);
         beh_vars{sub,1}=Behavior.FullOutput;
         labels{sub,1}=Behavior.FullOutputLabels;
         try
             parfor_progress;
         end
    end
    parfor_progress(0);
    toc
    subsample_beh_vars=array2table(squeeze(cell2nDMAT(beh_vars))','VariableNames',labels{1,1});
    save(SaveName,'subsample_beh_vars'); 
end
end 

function [ inFilesBySS,inFilesByRun ] = gradCPT_LoadBehavior(loadNames,SubSamples)
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
        inFilesBySS{n,1}.VTC=[inFilesBySS{n,1}.VTC;inFilesByRun{1,count}.VTC(SubSamples{count,1}==1,:)];
        inFilesBySS{n,1}.Zone=[inFilesBySS{n,1}.Zone;inFilesByRun{1,count}.Zone(SubSamples{count,1}==1,:)];
        inFilesBySS{n,1}.ZonePrime=[inFilesBySS{n,1}.ZonePrime;inFilesByRun{1,count}.ZonePrime(SubSamples{count,1}==1,:)];
        inFilesBySS{n,1}.VTCderiv=[inFilesBySS{n,1}.VTCderiv;inFilesByRun{1,count}.VTCderiv(SubSamples{count,1}==1,:)];        
        inFilesBySS{n,1}.ttt=[inFilesBySS{n,1}.ttt;inFilesByRun{1,count}.ttt(SubSamples{count,1}==1,:)];
        inFilesBySS{n,1}.response=[inFilesBySS{n,1}.response;inFilesByRun{1,count}.response(SubSamples{count,1}==1,:)];
        try
            inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;inFilesByRun{1,count}.data(SubSamples{count,1}==1,:)];
        catch
            if size(inFilesBySS{n,1}.data,2)<size(inFilesByRun{1,count}.data,2)
                inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;inFilesByRun{1,count}.data(SubSamples{count,1}==1,1:size(inFilesBySS{n,1}.data,2))];
            else
                temp=zeros(size(inFilesByRun{1,count}.data,1),size(inFilesBySS{n,1}.data,2));
                temp(:,1:size(inFilesByRun{1,count}.data,2))=inFilesByRun{1,count}.data;
                inFilesBySS{n,1}.data=[inFilesBySS{n,1}.data;temp(SubSamples{count,1}==1,:)];
            end    
        end    
        inFilesBySS{n,1}.bordertracker=[inFilesBySS{n,1}.bordertracker;inFilesByRun{1,count}.bordertracker(SubSamples{count,1}==1,:)];  
        inFilesBySS{n,1}.ComputeReward=inFilesByRun{1,count}.ComputeReward;
        inFilesBySS{n,1}.numTrials=inFilesBySS{n,1}.numTrials+inFilesByRun{1,count}.numTrials;
        count=count+1;
    end
end
end

function NewData=quikCellFormat(OldData)
    NewData=cell(1,size(OldData,2));
    for i = 1:size(OldData,2)
        NewData{1,i}=OldData(:,i);
    end
end