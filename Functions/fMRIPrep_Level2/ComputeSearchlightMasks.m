function [AnalysisParameters] = ComputeSearchlightMasks(fmriprep_table,ExperimentsDir,varargin)
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
[AnalysisName] = VariableSetter('AnalysisName',[''],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[slShape] = VariableSetter('slShape',[],varargin);
[slRadius] = VariableSetter('slRadius',[],varargin);
[slVoxelThreshold] = VariableSetter('slVoxelThreshold',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);


MaxSphereSLSize=[1,7;2,33;3,123;4,257;5,515;6,925;7,1419;8,2109;9,3017;10,4169;11,5575;12,7153;13,9171;14,11513;15,14147];
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SearchlightInfo';
%Allows you to set name for this particular analysis


%% Compile filepaths for input files for the analysis
try
    [filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName','raw3mm','TitleTextName','Select fMRI input for GLM:');
    filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
catch
    [filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','ActivationPatterns','TitleTextName','Select fMRI input for GLM:');
    filePaths_BaseTC=filePaths_BaseTC.WholeBrain;
end

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run','Session','Group'},'Perform analysis by subject or by run:',SingleSelect);
end
bySes=0;
bySS=0;
byGroup=0;
byRun=0;
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=find(fmriprep_table.numRuns_bySub)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
elseif strcmpi(SubjectOrRun,'Group')
    byGroup=1;
    useIndicies=find(fmriprep_table.numRuns_byGroup)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_byGroup;    
elseif strcmpi(SubjectOrRun,'Session')
    bySes=1;
    useIndicies=find(fmriprep_table.numRuns_bySes)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;      
else    
    byRun=1;
    useIndicies=[1:TotalRuns];
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;%% Apply graymatter mask
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
        gmAffix=['_gm',num2str4filename(gmMask,2)];
    end
end

if ~isempty(gmMask)
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end

%% Set searchligh setting
if isempty(slShape)
    SingleSelect=1; %Allows only a single value to be selected.
    [slShape] = uiNameSelect({'Sphere','Cube'},'Select searchlight shape:',SingleSelect);
end
if isempty(slRadius)
    slRadiusName=uiEnterName('3','Enter searchlight radius');
    slRadius=str2num(slRadiusName);
end

if strcmpi(slShape,'Sphere')
    maxVoxels=MaxSphereSLSize(slRadius,2);
else
    maxVoxels=(slRadius*2+1)^3;
end
if isempty(slVoxelThreshold)
    slVoxelThreshold=uiEnterName(num2str(round(maxVoxels*0.75)),['Enter searchlight radius ',newline,'Max # voxels: ',num2str(maxVoxels)]);
end

if isempty(AnalysisName)
    AnalysisName=uiEnterName(['SearchlightInfo_Shape-',slShape,'_Rad-',slRadiusName,'_Thresh-',slVoxelThreshold,gmAffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
end
slVoxelThreshold=str2num(slVoxelThreshold);
AnalysisParameters.slShape=slShape;
AnalysisParameters.slRadius=slRadius;
AnalysisParameters.slVoxelThreshold=slVoxelThreshold;
AnalysisParameters.SubjectOrRun=SubjectOrRun;
AnalysisParameters.gmMask=gmMask;

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
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end
    descript1='desc-SLinfo'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(1);
    SavePrefix=[ExperimentsDir,SaveDir,'/'];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    else  
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
            disp(['Skipping-- all files exist: ',SaveName]);
            continue
        end
    end
    %% Determine number of runs
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    %initialize variable names
    BaseTC=cell(1);
    brain_mask=cell(1);
    Use_gmMask=[];    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
        if ~isempty(fmriprep_table.Error{loadInd,1})
            continue
        end    
        %% set load paths and variable names
        %pull load paths
        try
            LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
        catch
            LoadPath_BaseTC=[];
        end     
        %Pull base timecourse data. If it doesn't exist, skip run.
        try
            variableInfo = who('-file', LoadPath_BaseTC);
            if any(ismember(variableInfo,'brain_mask'))                
                TempLoadData = load(LoadPath_BaseTC,'brain_mask'); 
            end
            brain_mask{count,1}=TempLoadData.brain_mask;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
            continue
        end
        if ~isempty(gmMask) && isempty(Use_gmMask)
            try
                TempLoadData = load(filePaths_gmMask{loadInd,1},'WM_probseg');
                Use_wmMask=single(TempLoadData.WM_probseg < 0.99);
            catch
                disp(['White matter mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                Use_wmMask=(brain_mask{count,1}*0)+1;
            end

            try
                TempLoadData = load(filePaths_gmMask{loadInd,1},'CSF_probseg');
                Use_csfMask=single(TempLoadData.CSF_probseg < 0.25);
            catch
                disp(['CSF mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                Use_csfMask=(brain_mask{count,1}*0)+1;
            end            
            try
                LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
                TempLoadData = load(LoadPath_gmMask,'GM_probseg');
                Use_gmMask=single(TempLoadData.GM_probseg > gmMask);
            catch
                disp(['Gray matter mask error-- mask not applied',LoadPath_gmMask]);
                Use_gmMask=(brain_mask{count,1}*0)+1;
            end
        elseif isempty(Use_gmMask)
            Use_gmMask=(brain_mask{count,1}*0)+1;
            Use_wmMask=(brain_mask{count,1}*0)+1;
            Use_csfMask=(brain_mask{count,1}*0)+1;
        end          
        count=count+1;
    end     
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    brain_mask=single(sum(single(cell2nDMAT(brain_mask)~=0),4)==size(brain_mask,1)).*Use_gmMask.*Use_csfMask.*Use_wmMask;
    [SLinds,SLnumVoxels,SLcoords] = Mask2Searchlight(brain_mask,slRadius,slVoxelThreshold,slShape);
    
    save(SaveNames{1,1},'SLinds','SLnumVoxels','SLcoords','AnalysisParameters');  
    
end