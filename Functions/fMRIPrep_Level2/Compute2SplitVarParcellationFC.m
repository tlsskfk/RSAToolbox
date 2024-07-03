function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = Compute2SplitVarParcellationFC(fmriprep_table,ExperimentsDir,varargin)
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
[SkipFC] = VariableSetter('SkipFC',0,varargin);
% Set Parcellation to run analysis on
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
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
[SplitAcrossRun] = VariableSetter('SplitAcrossRun',[],varargin);
[DefaultName] = VariableSetter('DefaultName',0,varargin);
[SplitVar] = VariableSetter('SplitVar',[],varargin);
[SplitLag] = VariableSetter('SplitLag',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
[NoParamEdit] = VariableSetter('NoParamEdit',0,varargin);
if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_SplitTC,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','SplitVarParcelFC','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_SplitTC.Properties.VariableNames;
        catch
            filePaths_SplitTC=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_SplitTC)
            for j = 1:height(filePaths_SplitTC)
                try
                    load(filePaths_SplitTC.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
                catch
                    continue
                end
                if ~isempty(AnalysisParameters)
                    break
                end
            end
        end
    end
end

if ~isempty(AnalysisParameters)
    ParamNames=fieldnames(AnalysisParameters);
    for i = 1:length(ParamNames)
        if isempty(AnalysisParameters.(ParamNames{i,1}))
            AnalysisParameters.(ParamNames{i,1})='NA_EMPTY';
        end
    end    
    if NoParamEdit==0
        SingleSelect=0; %Allows only a single value to be selected.
        [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect);    
        if ~isempty(EditParams)
            for i = 1:length(EditParams)
                AnalysisParameters.(EditParams{i,1})=[];
            end
        end
    end
    
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    gmMask=AnalysisParameters.gmMask;
    ParcelNames=AnalysisParameters.ParcelNames;
    SplitAcrossRun=AnalysisParameters.SplitAcrossRun;
    %NumReps=AnalysisParameters.NumReps;
    SplitVar=AnalysisParameters.SplitVar;
    SplitLag=AnalysisParameters.SplitLag;
    BaseTCName=AnalysisParameters.BaseTCName;
    RemoveStartVols=AnalysisParameters.RemoveStartVols;
    NormWithinRun=AnalysisParameters.NormWithinRun;
    BandPassFilter=AnalysisParameters.BandPassFilter;
    LoPassFilter=AnalysisParameters.LoPassFilter;
    HiPassFilter=AnalysisParameters.HiPassFilter;
    try
        fmriprep_table_name=AnalysisParameters.fmriprep_table_name;
    catch
        fmriprep_table_name=[];
    end
else
    AnalysisParameters=struct;
    Filters2Run=[];
    BaseTCName=[];
end
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
AnalysisParameters.fmriprep_table_name=fmriprep_table_name;
%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;
FilterAffix=['_filt'];

if ~strcmpi(BandPassFilter,'NA_EMPTY') || ~strcmpi(LoPassFilter,'NA_EMPTY') || ~strcmpi(HiPassFilter,'NA_EMPTY')
    if isempty(BandPassFilter) && isempty(LoPassFilter) && isempty(HiPassFilter)
        SingleSelect=0;
        [Filters2Run] = uiNameSelect({'None','HiPass','LoPass','BandPass'},'Select filters to apply:',SingleSelect);
    else
        Filters2Run=[];
        if ~isempty(HiPassFilter)
            Filters2Run=[Filters2Run;{'HiPass'}];
        end
        if ~isempty(LoPassFilter)
            Filters2Run=[Filters2Run;{'LoPass'}];
        end
        if ~isempty(BandPassFilter)
            Filters2Run=[Filters2Run;{'BandPass'}];
        end    
    end
else
    BandPassFilter=[];
    LoPassFilter=[];
    HiPassFilter=[];
    Filters2Run={'None'};
end
    
if ~strcmpi(Filters2Run,'None')
    if any(ismember(Filters2Run,'HiPass'))
        if isempty(HiPassFilter)
            HiPassFilter=uiEnterName('',['Enter hi pass filter threshold (in hz)']);
            HiPassFilter=str2num(HiPassFilter);
        end
        FilterAffix=[FilterAffix,'hp',num2str4filename(HiPassFilter)];
    end    
    if any(ismember(Filters2Run,'LoPass'))
        if isempty(LoPassFilter)
            LoPassFilter=uiEnterName('',['Enter lo pass filter threshold (in hz)']);
            LoPassFilter=str2num(LoPassFilter);
        end
        FilterAffix=[FilterAffix,'lp',num2str4filename(LoPassFilter)];
    end         
    if any(ismember(Filters2Run,'BandPass'))
        if isempty(BandPassFilter)
            tempBandPass=uiEnterName('',['Enter bandpass filter lower bound (in hz)']);
            BandPassFilter(1,1)=str2num(tempBandPass);
            tempBandPass=uiEnterName('',['Enter bandpass filter higher bound (in hz)']);
            BandPassFilter(1,2)=str2num(tempBandPass);
        end
        FilterAffix=[FilterAffix,'bp',num2str4filename(BandPassFilter(1,1)),'to',num2str4filename(BandPassFilter(1,2))];
    end 
end 
AnalysisParameters.HiPassFilter=HiPassFilter;
AnalysisParameters.LoPassFilter=LoPassFilter;
AnalysisParameters.BandPassFilter=BandPassFilter;
if strcmpi(FilterAffix,'_filt')
    FilterAffix='';
end
NormAffix=['_norm-'];
if isempty(NormWithinRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [NormWithinRun] = uiNameSelect({'zscore','none'},'Normalize time course data within run?',SingleSelect);
end
NormAffix=[NormAffix,NormWithinRun];
if strcmpi(NormWithinRun,'none')
    NormWithinRun='none';
end
AnalysisParameters.NormWithinRun=NormWithinRun;
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
    SplitAcrossRun=0;
end
if bySS==1
    if isempty(SplitAcrossRun)
        SingleSelect=1; %Allows only a single value to be selected.
        [SplitAcrossRun] = uiNameSelect({'Yes','No'},'Split across runs?',SingleSelect);
        if strcmpi(SplitAcrossRun,'Yes')
            SplitAcrossRun=1;
        else
            SplitAcrossRun=0;
        end
    end
end
AnalysisParameters.SplitAcrossRun=SplitAcrossRun;
RemoveVolAffix='_rm';
if isempty(RemoveStartVols)
    RemoveStartVols=uiEnterName('',['Remove start of each run? Leave blank for no.',newline,'Enter duration to remove in secs.']);
    RemoveStartVols=str2num(RemoveStartVols);
end
if strcmpi(RemoveStartVols,'NA_EMPTY')
    RemoveStartVols=[];
end
if isempty(RemoveStartVols)
    RemoveStartVols=[];
    RemoveVolAffix=[RemoveVolAffix,'0'];
else
    RemoveVolAffix=[RemoveVolAffix,num2str4filename(RemoveStartVols),'s'];
end
AnalysisParameters.RemoveStartVols=RemoveStartVols;

LagAffix='_lag-';
if isempty(SplitLag)
    SplitLag=uiEnterName('',['Enter split lag in secs.']);
    SplitLag=str2num(SplitLag);
end
if strcmpi(SplitLag,'NA_EMPTY')
    SplitLag=0;
end

LagAffix=[LagAffix,num2str4filename(SplitLag),'s'];

AnalysisParameters.SplitLag=SplitLag;

%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName',BaseTCName,'TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
AnalysisParameters.BaseTCName=BaseTCName;
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
[filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
filePaths_beh_events=filePaths_beh_events.(beh_eventsName);
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
        gmAffix=['_gm',num2str4filename(gmMask,2)];
    end
elseif isnumeric(gmMask)    
    gmAffix=['_gm',num2str4filename(gmMask,2)];
elseif strcmpi(gmMask,'NA_EMPTY')
    gmMask=[];
end
AnalysisParameters.gmMask=gmMask;
if ~isempty(gmMask)
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end
ConfoundNames=[];
for aNum=1:size(filePaths_ConfoundTCs,1)
    if ~isempty(filePaths_ConfoundTCs{aNum,1})            
        load(filePaths_ConfoundTCs{aNum,1},'ConfoundTCs');            
        ConfoundNames=unique([ConfoundNames;ConfoundTCs.Properties.VariableNames(:)]);
    end    
end 
EventNamesAll=[];
for aNum=1:size(filePaths_beh_events,1)
    if ~isempty(filePaths_beh_events{aNum,1})            
        load(filePaths_beh_events{aNum,1},'beh_events');            
        EventNamesAll=unique([EventNamesAll;beh_events.Properties.VariableNames(:)]);
    end    
end
if isempty(SplitVar)
    SplitVar=uiNameSelect([EventNamesAll(:);ConfoundNames(:)],'Select Split Var',1);
end
if iscell(SplitVar)
    SplitVar=SplitVar{1,1};
end
if any(ismember(ConfoundNames,SplitVar))
    SplitVarType='Confound';
elseif any(ismember(EventNamesAll,SplitVar))
    SplitVarType='Event';
end
AnalysisParameters.SplitVar=SplitVar;
SplitAffix=['_SplitBy',SplitVar];
if SplitAcrossRun == 1
    SplitAffix=[SplitAffix,'xRun'];
end

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='SplitVarParcelFC';

%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    tempName=[strrep(BaseTCName,'ResidTC_',''),SplitAffix,FilterAffix,NormAffix,RemoveVolAffix,LagAffix,gmAffix];
    if DefaultName~=1
        AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
    else
        AnalysisName=tempName;
    end    
end
 AnalysisNameTC=['ParcelTC_',strrep(BaseTCName,'ResidTC_',''),FilterAffix,NormAffix,RemoveVolAffix,gmAffix];
%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators

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
        [tempParcelFileNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','CoordParcels','AnalysisName',ParcelNames{i,1});
        AllParcelTable=[AllParcelTable,tempParcelFileNames];
    end
    ParcelNames=AllParcelTable.Properties.VariableNames(:);
elseif strcmpi(ParcelType,'IndvParcel')
    [AllParcelTable] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','IndvParcels');
    ParcelNames=uiNameSelect([AllParcelTable.Properties.VariableNames],'Select parcellations to run.');
    AllParcelTable=AllParcelTable(:,ParcelNames);
end
AnalysisParameters.ParcelNames=ParcelNames;
numParcels=length(ParcelNames);
%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];

if sum(ismember(ParcelNames,'WholeBrain'))==0 && strcmpi(ParcelType,'GroupParcel') 
    for i = 1:numParcels
        if i ==1            
            analysisMask=single(Parcels{i,1}.UseMask>0);
        else
            analysisMask=analysisMask+single(Parcels{i,1}.UseMask>0);
        end
        analysisMask=single(analysisMask>0);
    end
end

AnalysisParameters.RunDate=genDateString;

BaseAnalysisMask=analysisMask;
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    analysisMask=BaseAnalysisMask;
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    %% Set save directory and save name
    SaveDirTC=[fmriprep_table.matDir{dataInd,1},'/ParcelTC/',AnalysisNameTC,'/'];
    SaveNameTC = fmriprep_table.preproc_bold{dataInd,1};
    SaveNameTC=strrep(SaveNameTC,'.nii','');
    SaveNameTC=strrep(SaveNameTC,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end
    descript1='desc-parcelTCs'; %set file description name
    SaveNameTC=strrep(SaveNameTC,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNamesTC=cell(numParcels,1);
    RunParcel=ones(numParcels,1);     
    for parcelNum=1:numParcels
        SavePrefixTC=[ExperimentsDir,SaveDirTC,ParcelNames{parcelNum,1},'/'];    
        if ~exist(SavePrefixTC,'file')
            mkdir(SavePrefixTC);
            SaveNamesTC{parcelNum,1}=[SavePrefixTC,SaveNameTC,'.mat'];
        else  
            SaveNamesTC{parcelNum,1}=[SavePrefixTC,SaveNameTC,'.mat'];
            if exist(SaveNamesTC{parcelNum,1},'file')~=0 && Overwrite==0
                RunParcel(parcelNum,1)=0; 
            end
        end
    end
    if sum(RunParcel(:))==0
        disp(['Skipping ParcelTC-- all files exist: ',SaveNameTC]);
        SkipParcelTC=1;
    else
        SkipParcelTC=0;
    end
    
    SaveDir=[fmriprep_table.matDir{dataInd,1},'/',AnalysisType,'/',AnalysisName,'/'];
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
        SaveName=strrep(SaveName,'_run-1','');
    end
    descript1='desc-SplitVarFC'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold','');
    
    %% If analysis involves parcellations:
    SaveNames=cell(numParcels,1);
    for parcelNum=1:numParcels
        SavePrefix=[ExperimentsDir,SaveDir,ParcelNames{parcelNum,1},'/'];    
        if ~exist(SavePrefix,'file')
            mkdir(SavePrefix);
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
        else  
            SaveNames{parcelNum,1}=[SavePrefix,SaveName,'.mat'];
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
    Use_SplitVar=cell(1);
    Use_TrialNum=cell(1);
    Use_RunDur=cell(1);
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
        if SkipParcelTC==1
            try
                LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
                TempLoadData = load(LoadPath_ConfoundTCs,'ConfoundTCs');
                Use_RunDur{count,1}=height(TempLoadData.ConfoundTCs);
            catch
                disp('Error')
                continue
            end
        end
        if strcmpi(SplitVarType,'Confound')
            try
                LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
                TempLoadData = load(LoadPath_ConfoundTCs,'ConfoundTCs');
                Use_SplitVar{count,1}=TempLoadData.ConfoundTCs.(SplitVar);
                Use_TrialNum{count,1}=[1:length(Use_SplitVar{count,1})]';
            catch
                disp('Error')
                continue
            end           
        elseif strcmpi(SplitVarType,'Event')
            try               
                LoadPath_Events=filePaths_beh_events{loadInd,1};
                TempLoadData = load(LoadPath_Events,'beh_events');
                Use_SplitVar{count,1}=TempLoadData.beh_events.(SplitVar);
                Use_TrialNum{count,1}=TempLoadData.beh_events.TrialNum; 
                
            catch
                disp('Error')
                continue
            end
        end
        %Pull base timecourse data. If it doesn't exist, skip run.       
        if SkipParcelTC==0       
            try
                variableInfo = who('-file', LoadPath_BaseTC);
                if any(ismember(variableInfo,'boldTCs'))                
                    TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                    BaseTC{count,1}=TempLoadData.boldTCs; 
                elseif any(ismember(variableInfo,'Resids')) 
                    TempLoadData = load(LoadPath_BaseTC,'Resids','brain_mask');
                    BaseTC{count,1}=TempLoadData.Resids; 
                end
                Use_RunDur{count,1}=size(BaseTC{count,1},1);
                brain_mask{count,1}=TempLoadData.brain_mask;
            catch
                disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
                continue
            end   
            sampFreq=[];
            if ~isempty(LoPassFilter) || ~isempty(HiPassFilter) || ~isempty(BandPassFilter)
                try
                    sampFreq=1/fmriprep_table.TR(loadInd,1);
                catch
                    disp('TR unknown! Sample frequency for temporal filtering cannot be computed.'); 
                end
            end
            if ~isempty(sampFreq) && ~isempty(LoPassFilter)
                try
                    BaseTC{count,1}=lowpass(BaseTC{count,1},LoPassFilter,sampFreq);
                catch
                    disp('Lowpass filter failed!'); 
                end
            end    
            if ~isempty(sampFreq) && ~isempty(HiPassFilter)
                try
                    BaseTC{count,1}=highpass(BaseTC{count,1},HiPassFilter,sampFreq);
                catch
                    disp('Highpass filter failed!'); 
                end
            end   
            if ~isempty(sampFreq) && ~isempty(BandPassFilter)
                try
                    BaseTC{count,1}=bandpass(BaseTC{count,1},BandPassFilter,sampFreq);
                catch
                    disp('Bandpass filter failed!'); 
                end
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
        end
        count=count+1;
    end     
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    if SkipParcelTC==0 
        if isempty(analysisMask)
            analysisMask=Use_gmMask;
        else
            analysisMask=analysisMask.*Use_gmMask;
        end
        %% Run analysis here!!
        if ~isempty(RemoveStartVols)
            try
                tempTR=fmriprep_table.TR(dataInd,1);
                if iscell(tempTR)
                    tempTR=tempTR{1,1};
                end
                numVols2Remove=ceil(RemoveStartVols/tempTR);
            catch    
                numVols2Remove=[];
                disp('TR unknown! no vols removed at start of run');
            end
        else
            numVols2Remove=[];
        end
        [AllboldTCs,commonMask] = GLM_boldTCPrep(BaseTC,brain_mask,'NormWithinRun',NormWithinRun,'RemoveStartVols',numVols2Remove);
        brainSize=size(commonMask);
    end
    if ~isempty(SplitLag)
        try
            tempTR=fmriprep_table.TR(dataInd,1);
            if iscell(tempTR)
                tempTR=tempTR{1,1};
            end            
            UseSplitLag=ceil(SplitLag/tempTR);
        catch    
            UseSplitLag=0;
            disp('TR unknown! no SplitLag removed at start of run');
        end
    else
        UseSplitLag=0;
    end    

    % Resample SplitVarSize to volumes
	[Splits,SplitData] = ComputeVarSplits(Use_SplitVar,Use_TrialNum,SplitAcrossRun,Use_RunDur,UseSplitLag);
    if ~strcmpi(ParcelType,'GroupParcel')        
        tempParcelPaths=table2cell(AllParcelTable(dataInd,:));
        numParcels=size(tempParcelPaths,2);
        Parcels=cell(size(tempParcelPaths,2),1);
        for j = 1:size(tempParcelPaths,2)
            try
                Parcels{j,1}=load(tempParcelPaths{1,j},'UseMask','UseLabels'); 
            catch
                 Parcels{j,1}=[];
                 continue
            end
            if unique(Parcels{j,1}.UseMask(:)) == 0
                disp(['Skipping ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                continue
            end
        end
    end         
    for i = 1:numParcels
        if RunParcel(i,1)==1
            if ~strcmpi(ParcelType,'GroupParcel') 
                try
                    UseMask=Parcels{i,1}.UseMask;
                    UseLabels=Parcels{i,1}.UseLabels;  
                catch
                    continue
                end
            elseif strcmpi(ParcelNames{i,1},'WholeBrain')
                UseLabels=[];
                UseMask=analysisMask;
            else
                load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
                if ~ismember(brainSize,size(UseMask),'rows')
                    UseMask=imresize3(UseMask,brainSize,'nearest');
                end           
            end
            try
                [ParcelTC,ParcelVec,parcel_mask]=ComputeParcelTC(AllboldTCs,commonMask,UseMask,analysisMask,UseLabels);        
                save(SaveNamesTC{i,1},'ParcelTC','ParcelVec','parcel_mask','AnalysisParameters');      
            catch
                disp(['Error: ',ParcelNames{i,1}])
                continue
            end
        else
            load(SaveNamesTC{i,1},'ParcelTC');
        end
        UseLabels=ParcelTC.Properties.VariableNames(:);
        ParcelTC=table2array(ParcelTC);
        try
        [Split1_FCMat,Split1_FCVert] = ComputeFC(ParcelTC(Splits==1,:));  
        [Split2_FCMat,Split2_FCVert] = ComputeFC(ParcelTC(Splits==2,:));
        SplitDiff_FCMat=tanh(atanh(Split2_FCMat)-atanh(Split1_FCMat));
        SplitDiff_FCVert=mat2uppertriuvectormat(SplitDiff_FCMat);
        [Split1_FC,Split1_FCMat,Split1_FCMat_fZ,Split1_FCMat_ZNorm,Split1_FCMat_fZ_ZNorm,Split1_fZ_Degree,Split1_fZ_WholeCM,Split1_fZ_ZNorm_Degree,Split1_fZ_ZNorm_WholeCM] = FullParcelFC(Split1_FCMat,Split1_FCVert,UseLabels);
        [Split2_FC,Split2_FCMat,Split2_FCMat_fZ,Split2_FCMat_ZNorm,Split2_FCMat_fZ_ZNorm,Split2_fZ_Degree,Split2_fZ_WholeCM,Split2_fZ_ZNorm_Degree,Split2_fZ_ZNorm_WholeCM] = FullParcelFC(Split2_FCMat,Split2_FCVert,UseLabels);
        [SplitDiff_FC,SplitDiff_FCMat,SplitDiff_FCMat_fZ,SplitDiff_FCMat_ZNorm,SplitDiff_FCMat_fZ_ZNorm,SplitDiff_fZ_Degree,SplitDiff_fZ_WholeCM,SplitDiff_fZ_ZNorm_Degree,SplitDiff_fZ_ZNorm_WholeCM] = FullParcelFC(SplitDiff_FCMat,SplitDiff_FCVert,UseLabels);
        save(SaveNames{i,1},'Splits','SplitData','AnalysisParameters','Split1_FC','Split1_FCMat','Split1_FCMat_fZ','Split1_FCMat_ZNorm','Split1_FCMat_fZ_ZNorm','Split1_fZ_Degree','Split1_fZ_WholeCM','Split1_fZ_ZNorm_Degree','Split1_fZ_ZNorm_WholeCM',...
            'Split2_FC','Split2_FCMat','Split2_FCMat_fZ','Split2_FCMat_ZNorm','Split2_FCMat_fZ_ZNorm','Split2_fZ_Degree','Split2_fZ_WholeCM','Split2_fZ_ZNorm_Degree','Split2_fZ_ZNorm_WholeCM',...
            'SplitDiff_FC','SplitDiff_FCMat','SplitDiff_FCMat_fZ','SplitDiff_FCMat_ZNorm','SplitDiff_FCMat_fZ_ZNorm','SplitDiff_fZ_Degree','SplitDiff_fZ_WholeCM','SplitDiff_fZ_ZNorm_Degree','SplitDiff_fZ_ZNorm_WholeCM');
        catch
            continue
        end
    end        
    toc
end

%GroupAnalysis
for parcelNum=1:numParcels
    ParcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
        GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
    end
    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);     
    end
    if ~exist(GroupFiguresDir,'file')
        mkdir(GroupFiguresDir);     
    end
    if ~exist(GroupBrainMapsDir,'file')
        mkdir(GroupBrainMapsDir);     
    end     
    
    try
        load(['Parcellations/',ParcelName],'UseLabels','UseMask');
        UseLabels=UseLabels(:);
    catch
        try
            [compiledData] = CompileND(ExperimentsDir,fmriprep_table,'AnalysisType','CoordParcels','LoadVarName','ROInames','LoadVarFormat','Cell','DataType','Other','AnalysisName',ParcelName,'SubjectOrRun',SubjectOrRun,'TableOrArray','Table','VertOrMat','Matrix','ColumnCompile','All');
            UseLabels=compiledData(:,1,1);
            UseMask=[];
        catch
            load(['Parcellations/IndvParcels/',fmriprep_table_name,'/',ParcelName],'UseLabels','UseMask')
        end
    end
    [FCMats] = CompileND(ExperimentsDir,fmriprep_table,...
        'AnalysisType',AnalysisType,...
        'LoadVarName','Split1_FCMat_fZ',...
        'DataType','ByParcellation',...
        'LoadVarFormat','Table',...
        'AnalysisName',AnalysisName,...
        'TableOrArray','Array',...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelName);
    [Split1_Group_FC,Split1_Group_FCMat_Ts,Split1_Group_FCMat_Means,Split1_FC_Degree_Table] = GroupMatParcellationSummaryFigures(FCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,...
        'matName','Split1_FCz','SavePrefix','Split1f_');   
    [FCMats] = CompileND(ExperimentsDir,fmriprep_table,...
        'AnalysisType',AnalysisType,...
        'LoadVarName','Split2_FCMat_fZ',...
        'DataType','ByParcellation',...
        'LoadVarFormat','Table',...
        'AnalysisName',AnalysisName,...
        'TableOrArray','Array',...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelName);
    [Split2_Group_FC,Split2_Group_FCMat_Ts,Split2_Group_FCMat_Means,Split2_FC_Degree_Table] = GroupMatParcellationSummaryFigures(FCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,...
        'matName','Split2_FCz','SavePrefix','Split2_');   
     [FCMats] = CompileND(ExperimentsDir,fmriprep_table,...
        'AnalysisType',AnalysisType,...
        'LoadVarName','SplitDiff_FCMat_fZ',...
        'DataType','ByParcellation',...
        'LoadVarFormat','Table',...
        'AnalysisName',AnalysisName,...
        'TableOrArray','Array',...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelName);
    [SplitDiff_Group_FC,SplitDiff_Group_FCMat_Ts,SplitDiff_Group_FCMat_Means,SplitDiff_FC_Degree_Table] = GroupMatParcellationSummaryFigures(FCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,...
        'matName','SplitDiff_FCz','SavePrefix','SplitDiff_');      
    save([GroupAnalysisDir,'AnalysisParameters_FC'],'AnalysisParameters');
    save([GroupAnalysisDir,'Group_FC'],'Split1_Group_FC','Split2_Group_FC','SplitDiff_Group_FC');
    save([GroupAnalysisDir,'Group_FCMat'],'Split1_Group_FCMat_Ts','Split1_Group_FCMat_Means','Split1_FC_Degree_Table','Split2_Group_FCMat_Ts','Split2_Group_FCMat_Means','Split2_FC_Degree_Table','SplitDiff_Group_FCMat_Ts','SplitDiff_Group_FCMat_Means','SplitDiff_FC_Degree_Table');       
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

function [Splits,SplitData] = ComputeVarSplits(SplitVar,TrialNum,SplitAcrossRun,RunDurs,LagInVol)
    numRuns = size(SplitVar,1);
    VarByTrial=cell(numRuns,1);
    Splits=cell(numRuns,1);
    for i = 1:numRuns        
        tempSplitVar=SplitVar{i,1};
        Splits{i,1}=zeros(length(tempSplitVar),1);
        tempTrialNum=TrialNum{i,1};
        tempSplitVar(isnan(tempTrialNum),:)=[];
        tempTrialNum(isnan(tempTrialNum),:)=[];
        [~,UniqueInd]=unique(tempTrialNum);
        VarByTrial{i,1}=tempSplitVar(UniqueInd);
    end
    if SplitAcrossRun==1
        SplitVals=repmat(nanmedian(cell2mat(VarByTrial)),[numRuns,1]);
    else
        SplitVals=cellfun(@nanmedian,VarByTrial);
    end
    AllSplits=[];
    for i = 1:numRuns 
        Splits{i,1}(SplitVar{i,1}<SplitVals(i,1))=1;
        Splits{i,1}(SplitVar{i,1}>=SplitVals(i,1))=2;
        AllSplits=[AllSplits;Splits{i,1}];
        Splits{i,1}=imresize(Splits{i,1},[RunDurs{i,1},1],'nearest');
        if LagInVol > 0
            Splits{i,1}=[zeros(LagInVol,1);Splits{i,1}(1:end-LagInVol)];
        elseif LagInVol < 0
            Splits{i,1}=[Splits{i,1}(LagInVol+1:end);zeros(LagInVol,1);];
        end
    end
    SplitVar=cell2mat(SplitVar);
    SplitData=[nanmean(SplitVar(AllSplits==1),1),nanmean(SplitVar(AllSplits==2),1)];
    Splits=cell2mat(Splits);
end

function [FC,FCMat,FCMat_fZ,FCMat_ZNorm,FCMat_fZ_ZNorm,fZ_Degree,fZ_WholeCM,fZ_ZNorm_Degree,fZ_ZNorm_WholeCM] = FullParcelFC(FCMat,FCVert,UseLabels)
    
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

    [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(FCMat_fZ),'ThresholdType','rank','Threshold',0.5);
    rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
    for nameNum=1:length(rowNames50)
        rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
    end
    tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);

    [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(FCMat_fZ),'ThresholdType','rank','Threshold',0.25);
    rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
    for nameNum=1:length(rowNames25)
        rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
    end        
    tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);

    [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(FCMat_fZ),'ThresholdType','rank','Threshold',0.10);
    rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
    for nameNum=1:length(rowNames10)
        rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
    end        
    tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);

    [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(FCMat_fZ),'ThresholdType','rank','Threshold',0.05);
    rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
    for nameNum=1:length(rowNames05)
        rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
    end        
    tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);

    fZ_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
    fZ_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
    fZ_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];        

    [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(FCMat_fZ_ZNorm),'ThresholdType','rank','Threshold',0.5);
    rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
    for nameNum=1:length(rowNames50)
        rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
    end
    tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);

    [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(FCMat_fZ_ZNorm),'ThresholdType','rank','Threshold',0.25);
    rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
    for nameNum=1:length(rowNames25)
        rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
    end        
    tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);

    [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(FCMat_fZ_ZNorm),'ThresholdType','rank','Threshold',0.10);
    rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
    for nameNum=1:length(rowNames10)
        rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
    end        
    tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);

    [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(FCMat_fZ_ZNorm),'ThresholdType','rank','Threshold',0.05);
    rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
    for nameNum=1:length(rowNames05)
        rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
    end        
    tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);

    fZ_ZNorm_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
    fZ_ZNorm_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
    fZ_ZNorm_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];   
end       
        
            