function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeConfoundTC(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Template function for data processing from the fmriprep_table
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end
%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName','ConfoundTC',varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Run',varargin);


%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='confounds';

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
    descript1='desc-ConfoundTCs'; %set file description name
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
        LoadVars={'confounds_regressors'};
        try
            LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.funcDir{loadInd,1},fmriprep_table.(LoadVars{1,1}){dataInd,1}];
        catch
            LoadVars={'confounds_timeseries'};
            LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.funcDir{loadInd,1},fmriprep_table.(LoadVars{1,1}){dataInd,1}];
        end
        try
            [ ConfoundTCs ] = ConfoundLoad( LoadPaths{1,1},0 );
            fmriprep_table.numVol{loadInd,1}=size(ConfoundTCs,1);
            count=count+1;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,LoadPaths{1,1}]);
            if ~isempty(fmriprep_table.Error{loadInd,1})
                fmriprep_table.Error{loadInd,1}=[fmriprep_table.Error{loadInd,1},'_ConfoundLoad'];
            else
                fmriprep_table.Error{loadInd,1}='ConfoundLoad';
            end 
            continue
        end          
    end      
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ', SaveNames{1,1}]);
        continue
    end  
    
    %% correct for non steady state outliers
    regressor_names=ConfoundTCs.Properties.VariableNames;
    ConfoundTCs = fillmissing(ConfoundTCs,'linear');
    save(SaveNames{1,1},'ConfoundTCs','regressor_names');    
    toc
end
end 

function [ ConfoundRegressors ] = ConfoundLoad( FileName,combineNSSO )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ConfoundStruct=tdfread(FileName,'\t');
VarNames=fieldnames(ConfoundStruct);
NoJSON=0;
try
    jsonName=strrep(FileName,'.tsv','.json');
    jsonStruct = simple_JsonDecode(jsonName);
catch
    disp('No confound json file')
    NoJSON=1;
end
CSFCount=1;
WMCount=1;
CombinedCount=1;
for i=1:length(VarNames)
    VarName=VarNames{i,1};
    if NoJSON==0
        if isfield(jsonStruct,VarName)
            NewName=[];
            try
                tempField=jsonStruct.(VarName).Mask;
                if strcmpi(tempField,'CSF')
                    c=num2str(CSFCount);
                    padLength=4-length(c);
                    ZeroPad='';
                    if padLength > 0
                        for p = 1:padLength
                            ZeroPad=[ZeroPad,'0'];
                        end
                    end
                    NewName=['aCC_CSF_',ZeroPad,c];
                    CSFCount=CSFCount+1;
                elseif strcmpi(tempField,'WM')
                    c=num2str(WMCount);
                    padLength=4-length(c);
                    ZeroPad='';
                    if padLength > 0
                        for p = 1:padLength
                            ZeroPad=[ZeroPad,'0'];
                        end
                    end
                    NewName=['aCC_WM_',ZeroPad,c];
                    WMCount=WMCount+1;
                elseif strcmpi(tempField,'combined')
                    c=num2str(CombinedCount);
                    padLength=4-length(c);
                    ZeroPad='';
                    if padLength > 0
                        for p = 1:padLength
                            ZeroPad=[ZeroPad,'0'];
                        end
                    end
                    NewName=['aCC_combined_',ZeroPad,c];
                    CombinedCount=CombinedCount+1;
                end
                ConfoundStruct.(NewName)=ConfoundStruct.(VarName);
                ConfoundStruct=rmfield(ConfoundStruct,VarName);
                VarName=NewName;
                NewName=[];
            catch
                tempStruct=[];
            end
        end
    end
    if ~isnumeric(ConfoundStruct.(VarName))
        newMat=nan(size(ConfoundStruct.(VarName),1),1);
        removeInd=ismember(ConfoundStruct.(VarName)(:,1),'n');
        ConfoundStruct.(VarName)(removeInd,:)=[];
        tempMat=str2num(ConfoundStruct.(VarName));
        newMat(removeInd==0)=tempMat;
        ConfoundStruct.(VarName)=newMat;
    end
    if combineNSSO==1   
        if contains(VarName,'non_steady_state_outlier')
            if isfield(ConfoundStruct,'non_steady_state_outlier')
                ConfoundStruct.non_steady_state_outlier=ConfoundStruct.non_steady_state_outlier+ConfoundStruct.(VarName);
            else
                ConfoundStruct.non_steady_state_outlier=ConfoundStruct.(VarName);
            end
            ConfoundStruct=rmfield(ConfoundStruct,VarName);
        end
    end
end
ConfoundStruct.Linear=[1:length(ConfoundStruct.(VarNames{1,1}))]'/length(ConfoundStruct.(VarNames{1,1}));
ConfoundRegressors=struct2table(ConfoundStruct);

end