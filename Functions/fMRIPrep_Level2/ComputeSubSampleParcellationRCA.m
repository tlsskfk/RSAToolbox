function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeSubSampleParcellationRCA(fmriprep_table,ExperimentsDir,varargin)
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
[MatForm] = VariableSetter('MatForm',0,varargin);


if isempty(fmriprep_table_name)
    [~,fmriprep_table_names]=getFolderAndFileNames('fmriprep_table/');
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect(fmriprep_table_names,'Select fmriprep_table used:',SingleSelect);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
end
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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='subsample_RCA';
%Allows you to set name for this particular analysis

%% Compile filepaths for input files for the analysis
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','subsample_RSMs','TitleTextName',['Select subsample RSMs for',newline,'Representational Connectivity Analysis:']);
if isempty(AnalysisName)
    RCAName=strrep(RSMsName,'subsample_RSMs_',['subRCA_',fmriprep_table_name,'_']);
    RCAName=strrep(RSMsName,'subRSMs_',['subRCA_',fmriprep_table_name,'_']);
    AnalysisName=uiEnterName(RCAName,['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%% Identify parcellations availible and select ones to run analysis on.
ParcelNames=filePaths_RSMs.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_RSMs=filePaths_RSMs(:,ParcelNames);
numParcels=length(ParcelNames);

%% define save filepaths from the load filepaths 
[LoadDirs_RSMs,LoadNames_RSMs] = filepath2DirAndFileName(filePaths_RSMs);
oldDesc='desc-';
for i = 1:size(LoadNames_RSMs,1) %Identify old file description
    try
        [RSMsNameInfo] = fMRIPrep_Name2Info(LoadNames_RSMs{i,1}{1,1});
        oldDesc=[oldDesc,RSMsNameInfo.desc];
    catch
        continue
    end
    break
end
newDesc='desc-subRCA'; %set file description name
[SavePaths_RCA] = strrepCell(filePaths_RSMs,'/subsample_RSMs/',['/',AnalysisType,'/']);
[SavePaths_RCA] = strrepCell(SavePaths_RCA,['/',RSMsName,'/'],['/',AnalysisName,'/']);
[SavePaths_RCA] = strrepCell(SavePaths_RCA,['/',RSMsName,'/'],['/',AnalysisName,'/']);
[SaveDirs_RCA,SaveFileNames_RCA] = filepath2DirAndFileName(SavePaths_RCA);
[SaveFileNames_RCA] = strrepCell(SaveFileNames_RCA,oldDesc,newDesc);

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.RSMsName=RSMsName;
AnalysisParams.fmriprep_table_name=fmriprep_table_name;
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress


for parcelNum=1:numParcels
    tic
    parcelName=ParcelNames{parcelNum,1};
    disp(parcelName)
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
    end
    [RSMs,compiled_fmriprep_table,selectind_RSMs,UseLabels] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','RSMs',...
        'LoadVarFormat','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','subsample_RSMs',...
        'AnalysisName',RSMsName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',parcelName);
    if MatForm==1
        RSMsize=size(RSMs,1);
        RSMVertSize=(RSMsize^2-RSMsize)/2;
        numROIs=size(RSMs,3);
        numSubs=size(RSMs,4);
        numSS=size(RSMs,5);
        tempRSMs=RSMs;
        RSMs=zeros(RSMVertSize,numROIs,numSubs,numSS);
        for i = 1:numSubs
            for j = 1:numSS
                RSMs(:,:,i,j)=mat2uppertriuvectormat(tempRSMs(:,:,:,i,j));
            end
        end
        clearvars tempRSMs
    end
    [subsample_data_all,~,selectind_subsamples] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','subsample_data',...
        'LoadVarFormat','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','subsample_RSMs',...       
        'AnalysisName',RSMsName,...
        'SubjectOrRun',SubjectOrRun,...
        'TableOrArray','Array',...
        'ParcelName',parcelName);
    UseInd=find(selectind_RSMs==1 & selectind_subsamples==1);
    UseInd=UseInd(:);
    selectind_subsamples(selectind_RSMs==0)=[];
    RSMs(:,:,:,selectind_subsamples==0)=[];
    compiled_fmriprep_table(selectind_subsamples==0,:)=[];    
    if isempty(RSMs)
        continue
    end
    if ndims(RSMs)==5
        RSMs=permute(RSMs,[1,2,5,3,4]);
        numReps=size(RSMs,5);    
        subsample_data_all=permute(subsample_data_all,[1,3,2]);
    else
        RSMs=permute(RSMs,[1,2,4,3]);
        numReps=1;
    end
    numSubs=size(RSMs,4);
    numSS=size(RSMs,3);
    numROIs=size(RSMs,2); 

    AnalysisParams.UseTable=compiled_fmriprep_table;
    AnalysisParams.parcelName=parcelName;
    AnalysisParams.Parcellation=load(['Parcellations/',parcelName]);
    AnalysisParams.UseInd=UseInd;
    AnalysisParams.RSMPaths=strrepCell(filePaths_RSMs(UseInd,:).(parcelName),ExperimentsDir,'');
    AnalysisParams.RSMPaths=strrepCell(AnalysisParams.RSMPaths,'\','/');
    RSMs=real(RSMs); 
    for repNum= 1:numReps
        RFMats=zeros(numROIs,numSS,numSubs,'single');
        RCMats=zeros((numROIs^2-numROIs)/2,numSS,numSubs,'single');
        tempRSMs=RSMs(:,:,:,:,repNum);
        parfor subNum=1:numSubs
            [ OutVars ] = FastRCA( tempRSMs(:,:,:,subNum),'vertIn',1);
            RFMats(:,:,subNum)=single(OutVars.RF);
            RCMats(:,:,subNum)=single(OutVars.RCVert); 
        end    
        RFMats=permute(RFMats,[3,1,2]); %[Subsamples,ROIs,SS]
        RCMats=permute(RCMats,[3,1,2]); %[Subsamples,ROIs,SS]
        
        for loadInd = 1:length(UseInd)
            
            loadNum=UseInd(loadInd,1);
            subsample_data=subsample_data_all(:,loadInd);
            subRCA_RF=RFMats(:,:,loadInd);
            subRCA_RC=RCMats(:,:,loadInd);
            subRCA_RF_Corr = corr(subsample_data,subRCA_RF,'rows','pairwise');
            subRCA_RF_Corr = array2table(subRCA_RF_Corr,'VariableNames',UseLabels);
            
            subRCA_RC_Corr_Vert = corr(subsample_data,subRCA_RC,'rows','pairwise');
            [~,labelPairs2] = labels2uppertriuvectorlabels( UseLabels);
            [ subRCA_RC_Corr_Mat ] = vertRSM2SymRSM( subRCA_RC_Corr_Vert(:),length(UseLabels) );
%            [ ind,reverseInd ] = LabelPair2Ind( UseLabels,labelPairs2,'_2_' );
            
 %           subRCA_RC_Corr_Mat=coords2mat(cell2mat(ind),zeros(length(UseLabels),length(UseLabels)),subRCA_RC_Corr_Vert(:)) + coords2mat(cell2mat(reverseInd),zeros(length(UseLabels),length(UseLabels)),subRCA_RC_Corr_Vert(:));  
            subRCA_RC_Corr_Mat(eye(length(subRCA_RC_Corr_Mat))==1)=nan;
            subRCA_RC_Corr_Vert = array2table(subRCA_RC_Corr_Vert,'VariableNames',labelPairs2);

            if ~exist(SaveDirs_RCA.(parcelName){loadNum,1},'file')
                mkdir(SaveDirs_RCA.(parcelName){loadNum,1});
            end
            SaveName=[SaveDirs_RCA.(parcelName){loadNum,1},SaveFileNames_RCA.(parcelName){loadNum,1}];
            AnalysisParams.RCAPaths{loadInd,1}=SaveName;
            if repNum == 1
                save(SaveName,'subRCA_RF','subRCA_RC','subRCA_RF_Corr','subRCA_RC_Corr_Vert','subRCA_RC_Corr_Mat','AnalysisParams','subsample_data');
            else
                AppendSave(SaveName,'subRCA_RF',subRCA_RF,'Append',4);
                AppendSave(SaveName,'subsample_data',subsample_data,'Append',2);
                AppendSave(SaveName,'subRCA_RC',subRCA_RC,'Append',4);
                AppendSave(SaveName,'subRCA_RF_Corr',subRCA_RF_Corr,'Append',2);
                AppendSave(SaveName,'subRCA_RC_Corr_Vert',subRCA_RC_Corr_Vert,'Append',2);
                AppendSave(SaveName,'subRCA_RC_Corr_Mat',subRCA_RC_Corr_Mat,'Append',3);                
            end
        end
        
    end
    toc
end
end 