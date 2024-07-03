function [fmriprep_table] = ExtractAnatMasks(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Template function for data processing from the fmriprep_table

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Run',varargin);

Resample = VariableSetter( 'Resample',[],varargin);
VoxelSize = VariableSetter( 'VoxelSize',[],varargin);


%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='anat';

%Sets analysis name as boldTC_(VoxelSize)
if isempty(VoxelSize)
    SingleSelect=1;
    VoxelSize=uiNameSelect({'3mm','2mm','1mm'},'Select target voxel size (in mm^3).',SingleSelect);
end
AnalysisName=['probseg',VoxelSize];

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


%% BIDsTable loop: Iterate through fmriprep_table and perform analysis
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
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-probsegMasks'; %set file description name
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
        numRuns=fmriprep_table.numRuns(dataInd,1);
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
        LoadVars={'CSF_probseg','GM_probseg','WM_probseg'};
        LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.anatDir{loadInd,1},fmriprep_table.(LoadVars{1,1}){dataInd,1}];
        LoadPaths{1,2}=[ExperimentsDir,fmriprep_table.anatDir{loadInd,1},fmriprep_table.(LoadVars{1,2}){dataInd,1}];
        LoadPaths{1,3}=[ExperimentsDir,fmriprep_table.anatDir{loadInd,1},fmriprep_table.(LoadVars{1,3}){dataInd,1}];
        try
            CSF_probseg=load_nii(LoadPaths{1,1});
            CSF_probseg=single(CSF_probseg.img);
            maskDims=size(CSF_probseg);        
            count=count+1;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,LoadPaths{1,1}]);
            if ~isempty(fmriprep_table.Error{loadInd,1})
                fmriprep_table.Error{loadInd,1}=[fmriprep_table.Error{loadInd,1},'_CSFMaskLoad'];
            else
                fmriprep_table.Error{loadInd,1}='CSFMaskLoad';
            end 
            continue
        end       
        try
            GM_probseg=load_nii(LoadPaths{1,2});
            GM_probseg=single(GM_probseg.img);
            maskDims=size(GM_probseg);        
            count=count+1;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,2},newline,LoadPaths{1,2}]);
            if ~isempty(fmriprep_table.Error{loadInd,1})
                fmriprep_table.Error{loadInd,1}=[fmriprep_table.Error{loadInd,1},'_GMMaskLoad'];
            else
                fmriprep_table.Error{loadInd,1}='GMMaskLoad';
            end 
            continue
        end    
        try
            WM_probseg=load_nii(LoadPaths{1,3});
            WM_probseg=single(WM_probseg.img);
            maskDims=size(WM_probseg);        
            count=count+1;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,3},newline,LoadPaths{1,3}]);
            if ~isempty(fmriprep_table.Error{loadInd,1})
                fmriprep_table.Error{loadInd,1}=[fmriprep_table.Error{loadInd,1},'_WMMaskLoad'];
            else
                fmriprep_table.Error{loadInd,1}='WMMaskLoad';
            end 
            continue
        end            
    end  
    
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ',SaveNames{1,1}]);
        continue
    end  
    
    %% Resample brain mask (linear) to indicated size.
    if ~isempty(Resample)
        if ~ismember(maskDims,Resample)
            CSF_probseg=imresize3(CSF_probseg,Resample);
            GM_probseg=imresize3(GM_probseg,Resample);
            WM_probseg=imresize3(WM_probseg,Resample);
        end    
    end 
    save(SaveNames{1,1},'CSF_probseg','GM_probseg','WM_probseg');    
    toc
end
end 