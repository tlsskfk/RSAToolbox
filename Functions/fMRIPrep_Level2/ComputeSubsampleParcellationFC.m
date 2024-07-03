function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeSubsampleParcellationFC(fmriprep_table,ExperimentsDir,varargin)
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
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Subject or run level analysis. Will prompt request.
[Lag] = VariableSetter('Lag',[],varargin);

GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
AnalysisParams=struct;
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

if isempty(Lag)
    Lag=uiEnterName('',['Enter duration of TC lag in seconds',newline,'(usually 4-6sec)']);
    if isempty(LagName)
        Lag=0;
        LagName='Lag-0s';
    else
        Lag=str2num(Lag);
        LagName=['Lag-',num2str4filename(Lag),'s'];
    end
else
    LagName=['Lag-',num2str4filename(Lag),'s'];
end

[filePaths_SubSamples,~,SubSamplesName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','subsample','TitleTextName','Select sample sample type for RSM:');
filePaths_SubSamples=filePaths_SubSamples.(SubSamplesName);

%% Compile filepaths for input files for the analysis
[filePaths_ParcelTC,~,ParcelTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ParcelTC','TitleTextName','Select timecourses for functional connenctome:');
ParcelNames=filePaths_ParcelTC.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'WholeBrain'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_ParcelTC=filePaths_ParcelTC(:,ParcelNames);
numParcels=length(ParcelNames);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SubParcelFC';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    tempName=[AnalysisType,'_',LagName,'_',SubSamplesName,'_',ParcelTCName];
    tempName=strrep(tempName,'SubSample','sub');
    tempName=strrep(tempName,'subsample','sub');
    tempName=strrep(tempName,'_sub_','_');
    tempName=strrep(tempName,'__','_');
    AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.Parcellation=ParcelNames;
AnalysisParams.ParcelTCName=ParcelTCName;

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
    SaveDir=[fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-parcelFCs'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    
    for parcelNum=1:numParcels
        SaveNames=cell(1); 
        ParcelName=ParcelNames{parcelNum,1};
        SavePrefix=[ExperimentsDir,SaveDir,ParcelName,'/'];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
            SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        else  
            SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
            if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
                disp(['Skipping-- File exists: ',SaveNames{1,1}]);
                continue 
            end
        end
    
        %% Determine number of runs
        if bySS==1
            numRuns=fmriprep_table.numRuns(dataInd,1);
            for run=1:numRuns
                loadInd=dataInd+run-1;
                %skip previous errors
                if ~isempty(fmriprep_table.Error{loadInd,1})
                    continue
                end    
                %% set load paths and variable names
                %pull load paths
                LoadPath_subsamples{run,1}=filePaths_SubSamples{loadInd,1}; 
            end
        else
            numRuns=1;
        end 
        %% load input data 
        % If by Subject, iterate through runs and place data in cell
        if ~isempty(fmriprep_table.Error{dataInd,1})
            continue
        end    
        %% set load paths and variable names
        %pull load paths
        try               
            LoadPath_ParcelTC=filePaths_ParcelTC.(ParcelName){dataInd,1};
        catch
            LoadPath_ParcelTC=[];
        end
        %Pull base timecourse data. If it doesn't exist, skip run.
        if isempty(LoadPath_ParcelTC)
            continue
        end
        try
            TempLoadData = load(LoadPath_ParcelTC,'ParcelTC','AnalysisParams');
            ParcelTC=table2array(TempLoadData.ParcelTC); 
            AnalysisParams.ParcelTCParams=TempLoadData.AnalysisParams;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ParcelTC]);
            continue
        end               
        %% Run analysis here!!
        load(['Parcellations/',ParcelName],'UseLabels');
        [FCMat,FCVert] = ComputeFC(ParcelTC);  
        FCMat_fZ=array2table(atanh(FCMat),'VariableNames',UseLabels,'RowNames',UseLabels);
        FCMat_ZNorm=array2table(nan_zscore(FCMat,'pooled'),'VariableNames',UseLabels,'RowNames',UseLabels);
        FCMat_fZ_ZNorm=array2table(nan_zscore(atanh(FCMat),'pooled'),'VariableNames',UseLabels,'RowNames',UseLabels);
        FCMat=array2table(FCMat,'VariableNames',UseLabels,'RowNames',UseLabels);        
        FCVert_fZ=atanh(FCVert);
        FCVert_ZNorm=nan_zscore(FCVert,'pooled');
        FCVert_fZ_ZNorm=nan_zscore(FCVert_fZ,'pooled');        

        [labelPairs1,labelPairs2,labelPairsCell]=labels2uppertriuvectorlabels(UseLabels);
        [ FC_ind,FC_reverseInd ] = LabelPair2Ind( UseLabels,labelPairs1,'_2_' );
        FC=array2table([FCVert,FCVert_fZ,FCVert_ZNorm,FCVert_fZ_ZNorm],'VariableNames',{'FC','FC_fisherZ','FC_ZNorm','FC_fisherZ_ZNorm'},'RowNames',labelPairs2);
        FC=[FC,cell2table(labelPairsCell,'VariableNames',{'ROI1','ROI2'}),cell2table([labelPairs1,labelPairs2],'VariableNames',{'LabelPairs1','LabelPairs2'}),cell2table([ FC_ind,FC_reverseInd ],'VariableNames',{'MatInd','ReverseMatInd'})];                    
        save(SaveNames{1,1},'FC','FCMat','FCMat_fZ','FCMat_ZNorm','FCMat_fZ_ZNorm','AnalysisParams');
    end     
toc    
end
end 