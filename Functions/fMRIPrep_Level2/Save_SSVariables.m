function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = Save_SSVariables(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end
%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Subject',varargin);
[TRTSS] = VariableSetter('TRTSS',[],varargin);
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SSVars';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['SSVars_All'],['Enter name for ',AnalysisType,newline,'analysis below:']);
end

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
    
    SaveDir=[ExperimentsDir,fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-SS_Vars'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    if ~exist(SaveDir,'file')
        mkdir(SaveDir);
    end    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    loadName=[ExperimentsDir,fmriprep_table.experiment{dataInd,1},'/SS_Vars.mat'];
    try
        TempLoadData = load(loadName,'SS_Vars');
    catch
        disp('Skipping run-- input file or variable doesnt exist:');
        continue
    end  
    
    if TRTSS==1
        ssID=fmriprep_table.orig_IDs{dataInd,1};
        try
            SS_Vars=TempLoadData.SS_Vars(ssID,:);
        catch
            if length(ssID)==8
                ssID=[ssID,'_01'];
            else
                ssID=ssID(1,1:8);
            end
            try
                SS_Vars=TempLoadData.SS_Vars(ssID,:);
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ']);
                continue
            end
        end
        
  
    else        
        try
            SS_Vars=TempLoadData.SS_Vars(fmriprep_table.orig_IDs{dataInd,1},:);
        catch
            try
                SS_Vars=TempLoadData.SS_Vars(fmriprep_table.sub{dataInd,1},:);
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ']);
                continue
            end
        end
    end
    save([SaveDir,SaveName],'SS_Vars');
    toc
end
end 