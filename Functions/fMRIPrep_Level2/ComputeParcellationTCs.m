function [AnalysisParameters] = ComputeParcellationTCs(fmriprep_table,ExperimentsDir,varargin)
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
% Set Parcellation to run analysis on
[ParcelNames] = VariableSetter('Parcellation',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
% Remove N volumes at the beginning of run. Enter time in secs
[RemoveStartVols] = VariableSetter('RemoveStartVols',[],varargin);
% Z normalize within run
[NormWithinRun] = VariableSetter('NormWithinRun',[],varargin);
% BandPass filter
[BandPassFilter] = VariableSetter('BandPassFilter',[],varargin);
% LoPass filter
[LoPassFilter] = VariableSetter('LoPassFilter',[],varargin);
% HiPass filter
[HiPassFilter] = VariableSetter('HiPassFilter',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
AnalysisParameters=struct;
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
if isempty(fmriprep_table_name)
[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
FilterAffix=['_filt'];
if isempty(BandPassFilter) && isempty(LoPassFilter) && isempty(HiPassFilter)
    SingleSelect=0;
    [Filters2Run] = uiNameSelect({'None','HiPass','LoPass','BandPass'},'Select filters to apply:',SingleSelect);
else
    Filters2Run=[];
    FilterAffix=[FilterAffix,'None'];
end

if ~isempty(Filters2Run) && ~strcmpi(Filters2Run,'None')
    if any(ismember(Filters2Run,'HiPass'))
        HiPassFilter=uiEnterName('',['Enter hi pass filter threshold (in hz)']);
        HiPassFilter=str2num(HiPassFilter);
        FilterAffix=[FilterAffix,'hp',num2str4filename(HiPassFilter)];
    end    
    if any(ismember(Filters2Run,'LoPass'))
        LoPassFilter=uiEnterName('',['Enter lo pass filter threshold (in hz)']);
        LoPassFilter=str2num(LoPassFilter);
        FilterAffix=[FilterAffix,'lp',num2str4filename(LoPassFilter)];
    end         
    if any(ismember(Filters2Run,'BandPass'))
        tempBandPass=uiEnterName('',['Enter bandpass filter lower bound (in hz)']);
        BandPassFilter(1,1)=str2num(tempBandPass);
        tempBandPass=uiEnterName('',['Enter bandpass filter higher bound (in hz)']);
        BandPassFilter(1,2)=str2num(tempBandPass);
        FilterAffix=[FilterAffix,'bp',num2str4filename(BandPassFilter(1,1)),'to',num2str4filename(BandPassFilter(1,2))];
    end 
end    
NormAffix=['_norm-'];

if isempty(NormWithinRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [NormWithinRun] = uiNameSelect({'zscore','none'},'Normalize time course data within run?',SingleSelect);
    NormAffix=[NormAffix,NormWithinRun];
    if strcmpi(NormWithinRun,'none')
        NormWithinRun='';
    end
end

SubjectAffix='';
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    if length(unique(fmriprep_table.session))>1
        useIndicies=find(fmriprep_table.numRuns_bySes)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;        
    else
        useIndicies=find(fmriprep_table.numRuns_bySub)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    end
    SubjectAffix='_bySS';
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end

RemoveVolAffix='_rm';
if isempty(RemoveStartVols)
    RemoveStartVols=uiEnterName('',['Remove start of each run? Leave blank for no.',newline,'Enter duration to remove in secs.']);
    if isempty(RemoveStartVols)
        RemoveStartVols=[];
        RemoveVolAffix=[RemoveVolAffix,'0'];
    else
        RemoveStartVols=str2num(RemoveStartVols);
        RemoveVolAffix=[RemoveVolAffix,num2str4filename(RemoveStartVols),'s'];
    end
end

%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
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

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ParcelTC';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['ParcelTC_',strrep(BaseTCName,'ResidTC_',''),FilterAffix,NormAffix,RemoveVolAffix,gmAffix,SubjectAffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
end
%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators

if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

%% Parcellation processing
%If necessary, select parcellations to run analyses on. Output of this
%process are the variables:
    %ParcelNames: (N by 1 cell) containing the names of the parcellations 
    %numParcels: Number of parcellations (N) in ParcelNames. 
if ischar(ParcelNames) && ~strcmpi(ParcelNames,'NA_EMPTY')
    ParcelNames={ParcelNames};
end    
if isempty(ParcelType)
    ParcelType=uiNameSelect([{'ParcelCoords','GroupParcel','IndvParcel'}],'Select parcellations to run.');
end
AnalysisParameters.ParcelType=ParcelType; 
if strcmpi(ParcelType,'GroupParcel')
    if isempty(ParcelNames)
        ParcelNames=cell(1);
        ParcelDirInfo=dir('Parcellations/');
        count = 1;
        for i = 1:size(ParcelDirInfo,1)
            if ParcelDirInfo(i).isdir==0
                ParcelNames{count,1}=strrep(ParcelDirInfo(i).name,'.mat','');
                count=count+1;
            end
        end
        AllParcelNames=ParcelNames;
        ParcelNames=uiNameSelect([{'Searchlight'};ParcelNames],'Select parcellations to run.');
    elseif iscell(ParcelNames)
        ParcelNames=ParcelNames(:);
    end
    filePaths_Parcels=[];
elseif strcmpi(ParcelType,'ParcelCoords')
    [ParcelNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','CoordParcels','GetAnalysisNames',1);
    ParcelNames=uiNameSelect([ParcelNames],'Select parcellations to run.');
    AllParcelTable=[];
    for i = 1:length(ParcelNames)
        [tempParcelFileNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','CoordParcels','AnalysisName',ParcelNames{i,1});
        [tempParcelFileNames] = sub2runFilePath(fmriprep_table.numRuns_bySub,tempParcelFileNames);
        AllParcelTable=[AllParcelTable,tempParcelFileNames];
    end
    ParcelNames=AllParcelTable.Properties.VariableNames(:);
elseif strcmpi(ParcelType,'IndvParcel')
    [AllParcelTable] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','IndvParcels');
    ParcelNames=uiNameSelect([AllParcelTable.Properties.VariableNames],'Select parcellations to run.');
    AllParcelTable=AllParcelTable(:,ParcelNames);
end    
    
numParcels=length(ParcelNames);

%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];

if sum(ismember(ParcelNames,'WholeBrain'))==0
    for i = 1:numParcels
        parcelName=ParcelNames{i};
        try
            load(['Parcellations/',parcelName],'UseLabels','UseMask');
            UseLabels=UseLabels(:);
        catch
            try
                [compiledData] = CompileND(ExperimentsDir,fmriprep_table,'AnalysisType','CoordParcels','LoadVarName','ROInames','LoadVarFormat','Cell','DataType','Other','AnalysisName',parcelName,'SubjectOrRun',SubjectOrRun,'TableOrArray','Table','VertOrMat','Matrix','ColumnCompile','All');
                UseLabels=compiledData(:,1,1);
                UseMask=[];
            catch
                load(['Parcellations/IndvParcels/',fmriprep_table_name,'/',parcelName],'UseLabels','UseMask')
            end
        end
        if i ==1            
            analysisMask=single(UseMask>0);
        else
            analysisMask=analysisMask+single(UseMask>0);
        end
        analysisMask=single(analysisMask>0);
    end
end

BaseAnalysisMask=analysisMask;

if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end
AnalysisParameters.AnalysisType=AnalysisType;
AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.Parcellation=ParcelNames;
AnalysisParameters.gmMask=gmMask;
AnalysisParameters.RemoveStartVols=RemoveStartVols;
AnalysisParameters.NormWithinRun=NormWithinRun;
AnalysisParameters.BandPassFilter=BandPassFilter;
AnalysisParameters.LoPassFilter=LoPassFilter;
AnalysisParameters.HiPassFilter=HiPassFilter;
AnalysisParameters.BaseTCName=BaseTCName;


%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    try
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    analysisMask=BaseAnalysisMask;
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
    descript1='desc-parcelTCs'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(numParcels,1);
    RunParcel=ones(numParcels,1);     
    for parcelNum=1:numParcels
        SavePrefix=[ExperimentsDir,SaveDir,ParcelNames{parcelNum,1},'/'];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
        else  
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
            if exist(SaveNames{parcelNum,1},'file')~=0 && Overwrite==0
                RunParcel(parcelNum,1)=0; 
            end
        end
    end
    if sum(RunParcel(:))==0
        disp(['Skipping-- all files exist: ',SaveName]);
        continue
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
            if any(ismember(variableInfo,'boldTCs'))                
                TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                BaseTC{count,1}=TempLoadData.boldTCs; 
            elseif any(ismember(variableInfo,'Resids')) 
                TempLoadData = load(LoadPath_BaseTC,'Resids','brain_mask');
                BaseTC{count,1}=TempLoadData.Resids; 
            end
            brain_mask{count,1}=TempLoadData.brain_mask;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
            continue
        end        
        if ~isempty(gmMask) && isempty(Use_gmMask)
            try
                LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
            catch
                LoadPath_gmMask=[];
            end
            try
                TempLoadData = load(LoadPath_gmMask,'GM_probseg');
                Use_gmMask=single(TempLoadData.GM_probseg > gmMask);
            catch
                disp(['Gray matter mask error-- mask not applied',LoadPath_gmMask]);
                Use_gmMask=(brain_mask{count,1}*0)+1;
            end
        elseif isempty(Use_gmMask)
            Use_gmMask=(brain_mask{count,1}*0)+1;
        end       
        count=count+1;
    end     
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    if isempty(analysisMask)
        analysisMask=Use_gmMask;
    else
        try
            analysisMask=analysisMask.*Use_gmMask;
        catch
            analysisMask=imresize3(analysisMask,size(Use_gmMask),'nearest');
            analysisMask=analysisMask.*Use_gmMask;
        end
    end
    %% Run analysis here!!
    if ~isempty(RemoveStartVols)
        try
            if iscell(fmriprep_table.TR)
                numVols2Remove=ceil(RemoveStartVols/fmriprep_table.TR{dataInd,1});           
            else    
                numVols2Remove=ceil(RemoveStartVols/fmriprep_table.TR(dataInd,1));
            end
        catch    
            numVols2Remove=[];
            disp('TR unknown! no vols removed at start of run');
        end
    else
        numVols2Remove=[];
    end
    [AllboldTCs,commonMask] = GLM_boldTCPrep(BaseTC,brain_mask,'NormWithinRun',NormWithinRun,'RemoveStartVols',numVols2Remove);
    sampFreq=[];
    if ~isempty(LoPassFilter) || ~isempty(HiPassFilter) || ~isempty(BandPassFilter)
        try
            if any(iscell(fmriprep_table.TR))
                sampFreq=1/fmriprep_table.TR{dataInd,1};
            else    
                sampFreq=1/fmriprep_table.TR(dataInd,1);
            end
        catch
            disp('TR unknown! Sample frequency for temporal filtering cannot be computed.'); 
        end
    end
    if ~isempty(sampFreq) && ~isempty(LoPassFilter)
        try
            AllboldTCs=lowpass(AllboldTCs,LoPassFilter,sampFreq);
        catch
            disp('Lowpass filter failed!'); 
        end
    end    
    if ~isempty(sampFreq) && ~isempty(HiPassFilter)
        try
            AllboldTCs=highpass(AllboldTCs,HiPassFilter,sampFreq);
        catch
            disp('Highpass filter failed!'); 
        end
    end   
    if ~isempty(sampFreq) && ~isempty(BandPassFilter)
        try
            AllboldTCs=bandpass(AllboldTCs,BandPassFilter,sampFreq);
        catch
            disp('Bandpass filter failed!'); 
        end
    end       
    brainSize=size(commonMask);
    for i = 1:numParcels
        if strcmpi(ParcelNames{i,1},'WholeBrain')
            UseLabels=[];
            UseMask=analysisMask;
        elseif strcmpi(ParcelType,'GroupParcel') 
            Parcels=load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
            if ~ismember(brainSize,size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize,'nearest');
            end                        
        else
            Parcels=load(AllParcelTable.(ParcelNames{i,1}){dataInd,1},'UseMask','UseLabels');          
            if unique(Parcels.UseMask(:)) == 0
                disp(['Skipping ',AllParcelTable.(ParcelNames{i,1}){dataInd,1},' -- ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                continue
            end
        end
        try
            [ParcelTC,ParcelVec,parcel_mask]=ComputeParcelTC(AllboldTCs,commonMask,Parcels.UseMask,analysisMask,Parcels.UseLabels);        
            save(SaveNames{i,1},'ParcelTC','ParcelVec','parcel_mask','AnalysisParameters');      
        catch
            disp(['Error: ',ParcelNames{i,1}])
            continue
        end
    end        
    toc
    catch
        disp('error')
    end
end
end 

function [ParcelTC,ParcelVec,parcel_mask]=ComputeParcelTC(boldTCs,TCMask,ParcelMask,AnatMask,UseLabels)
if ~isempty(UseLabels)    
    ParcelTC=nan(size(boldTCs,1),length(UseLabels),'single');
end
parcel_mask=single(ParcelMask.*AnatMask.*single(TCMask~=0));
VertMask=parcel_mask(:);
VertMask(TCMask(:)==0,:)=[];
boldTCs(:,VertMask==0)=[];
ParcelVec=VertMask(VertMask~=0);
parcelNums=unique(ParcelVec)';
if ~isempty(UseLabels)  
    for i = parcelNums
        try
            ParcelTC(:,i)=mean(boldTCs(:,ParcelVec==i),2);    
        catch
            continue
        end
    end
    ParcelTC=array2table(ParcelTC,'VariableNames',UseLabels);
else
    ParcelTC=boldTCs;
end
end