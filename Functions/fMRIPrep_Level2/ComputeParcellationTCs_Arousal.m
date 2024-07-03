function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeParcellationTCs_Arousal(fmriprep_table,ExperimentsDir,varargin)
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
[SubjectOrRun] = VariableSetter('SubjectOrRun',['Run'],varargin);
% Set Parcellation to run analysis on
[Parcellation] = VariableSetter('Parcellation',['ForthVentrical_3mmMNI'],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
[roiSize] = VariableSetter('roiSize',[60],varargin);
% Remove N volumes at the beginning of run. Enter time in secs
[RemoveStartVols] = VariableSetter('RemoveStartVols',[0],varargin);
% Z normalize within run
[NormWithinRun] = VariableSetter('NormWithinRun',['none'],varargin);
% BandPass filter
[BandPassFilter] = VariableSetter('BandPassFilter',['none'],varargin);
% LoPass filter
[LoPassFilter] = VariableSetter('LoPassFilter',['none'],varargin);
% HiPass filter
[HiPassFilter] = VariableSetter('HiPassFilter',['none'],varargin);

AnalysisParams=struct;
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
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
elseif strcmpi(NormWithinRun,'none')
    NormAffix='';
    NormWithinRun='';
end

SubjectAffix='';
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    useIndicies=dataInd;
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
[filePaths_ConfoundTCs,~,ConfoundTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_ConfoundTCs=filePaths_ConfoundTCs.(ConfoundTCName);
[filePaths_beh_events,~,beh_eventsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh','AnalysisName','events');
filePaths_beh_events=filePaths_beh_events.(beh_eventsName);
[filePaths_behVars,~,behVarsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','beh_vars','AnalysisName','Overall');
filePaths_behVars=filePaths_behVars.(behVarsName);
gmAffix=[];
% if isempty(gmMask)
%     gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
%     if ~isempty(gmMask)
%         gmMask=str2num(gmMask);
%         gmAffix=['_gm',num2str4filename(gmMask,2)];
%     end
% elseif strcmpi(gmMask,'none')
%     gmMask=[];
% end

% if ~isempty(gmMask)
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
% end

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ArousalTC';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['ArousalTC_',strrep(BaseTCName,'ResidTC_',''),FilterAffix,NormAffix,RemoveVolAffix,gmAffix,SubjectAffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
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
if ischar(Parcellation)
    ParcelNames={Parcellation};
elseif ~iscell(Parcellation)
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
    ParcelNames=uiNameSelect([{'WholeBrain'};ParcelNames],'Select parcellations to run.');
else
    ParcelNames=Parcellation(:);
end
numParcels=length(ParcelNames);

%Create mask that excludes voxels that aren't included in any parcellation.
%Helpful to save time and minimize computation.
analysisMask=[];

if sum(ismember(ParcelNames,'WholeBrain'))==0
    for i = 1:numParcels
        load(['Parcellations/',ParcelNames{i,1}],'UseMask');
        if i ==1            
            analysisMask=single(UseMask>0);
        else
            analysisMask=analysisMask+single(UseMask>0);
        end
        analysisMask=single(analysisMask>0);
    end
end
BaseAnalysisMask=analysisMask;
if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.Parcellation=ParcelNames;
AnalysisParams.gmMask=gmMask;
AnalysisParams.RemoveStartVols=RemoveStartVols;
AnalysisParams.NormWithinRun=NormWithinRun;
AnalysisParams.BandPassFilter=BandPassFilter;
AnalysisParams.LoPassFilter=LoPassFilter;
AnalysisParams.HiPassFilter=HiPassFilter;
AnalysisParams.BaseTCName=BaseTCName;


%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    tic
    %% Display progress
    analysisMask=BaseAnalysisMask;
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
        SaveName=strrep(SaveName,'_run-1','');
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
        try
            LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
        catch
            LoadPath_gmMask=[];
        end
        try
            TempLoadData = load(LoadPath_gmMask,'CSF_probseg');
            tempCSFMask=TempLoadData.CSF_probseg.*analysisMask.*single(brain_mask{count,1}~=0);
            useProb=sort(tempCSFMask(:),'descend');
            useProb=useProb(roiSize,1);
            Use_gmMask=single(TempLoadData.CSF_probseg >= useProb);
            disp(['CSF thresh: ', num2str(useProb,2)]);
        catch
            disp(['Gray matter mask error-- mask not applied',LoadPath_gmMask]);
            Use_gmMask=(brain_mask{count,1}*0)+1;
        end
        LoadPath_Events=filePaths_beh_events{loadInd,1};
        LoadPath_behVars=filePaths_behVars{loadInd,1};
        TempLoadData=load(LoadPath_Events,'beh_events');
        beh_events=TempLoadData.beh_events;
        LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
        TempLoadData=load(LoadPath_ConfoundTCs,'ConfoundTCs');
        ConfoundTCs=TempLoadData.ConfoundTCs;
        globalSignal=ConfoundTCs.global_signal/nanmean(ConfoundTCs.global_signal,1)*100;
        wmSignal=ConfoundTCs.white_matter/nanmean(ConfoundTCs.white_matter,1)*100;
        beh_events.VigilenceTC_global=nan(height(beh_events),1);
        beh_events.VigilenceTC_CSF4=nan(height(beh_events),1);
        TempLoadData=load(LoadPath_behVars,'beh_vars');
        beh_vars=TempLoadData.beh_vars;
        count=count+1;
        
    end     
    runDur=fmriprep_table.fmriDur{dataInd,1};
    TR=fmriprep_table.TR(dataInd,1);
    maxTime=max(beh_events.TrialOnsetTime);
    if maxTime<=runDur
        endInd=size(beh_events.TrialOnsetTime,1);
    else
        [~,endInd]=min(abs(beh_events.TrialOnsetTime-runDur));
    end
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    if isempty(analysisMask)
        analysisMask=Use_gmMask;
    else
        analysisMask=analysisMask.*Use_gmMask;
    end
    %% Run analysis here!!
    if ~isempty(RemoveStartVols)
        try
            numVols2Remove=ceil(RemoveStartVols/fmriprep_table.TR(dataInd,1));
        catch    
            numVols2Remove=[];
            disp('TR unknown! no vols removed at start of run');
        end
    else
        numVols2Remove=[];
    end
    
    [AllboldTCs,commonMask] = GLM_boldTCPrep(BaseTC,brain_mask,'NormWithinRun',NormWithinRun,'RemoveStartVols',numVols2Remove);
    sampFreq=1/TR;  
    numVol=size(AllboldTCs,1);
    AllboldTCs=AllboldTCs./repmat(nanmean(AllboldTCs,1),[size(AllboldTCs,1),1])*100;
    brainSize=size(commonMask);
    for i = 1:numParcels
        if strcmpi(ParcelNames{i,1},'WholeBrain')
            UseLabels=[];
            UseMask=analysisMask;
        else
            load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
            if ~ismember(brainSize,size(UseMask),'rows')
                UseMask=imresize3(UseMask,brainSize,'nearest');
            end           
        end
        try
            [ParcelTC,~,parcelMask]=ComputeParcelTC(AllboldTCs,commonMask,UseMask,analysisMask,UseLabels);  
            %FourthTC=zscore(ParcelTC.FourthVentrical);
            disp(['FV size: ',num2str(sum(parcelMask(:))),' voxels']);
            FourthTC=ParcelTC.FourthVentrical;
            WindowSize=ceil(60/TR);
            ResultMatGlobal=nan(numVol-WindowSize+1,numVol);
            ResultMatGlobal2=nan(numVol-WindowSize+1,numVol);
            ResultMatCSF=nan(numVol-WindowSize+1,numVol);
            ResultMatCSFcontrol=ResultMatCSF;
            quarters=round(1:numVol/4:numVol);
            FVQuart1=FourthTC(1:quarters(1,2)-1,1);
            FVQuart2=FourthTC(quarters(1,2):quarters(1,3)-1,1);
            FVQuart3=FourthTC(quarters(1,3):quarters(1,4)-1,1);
            FVQuart4=FourthTC(quarters(1,4):end,1);
            
            GlobalQuart1=globalSignal(1:quarters(1,2)-1,1);
            GlobalQuart2=globalSignal(quarters(1,2):quarters(1,3)-1,1);
            GlobalQuart3=globalSignal(quarters(1,3):quarters(1,4)-1,1);
            GlobalQuart4=globalSignal(quarters(1,4):end,1);           
            wmQuart1=wmSignal(1:quarters(1,2)-1,1);
            wmQuart2=wmSignal(quarters(1,2):quarters(1,3)-1,1);
            wmQuart3=wmSignal(quarters(1,3):quarters(1,4)-1,1);
            wmQuart4=wmSignal(quarters(1,4):end,1);     
            
            BehResults.FV_Q1=bandpower(FVQuart1,sampFreq,[0.03,0.07]);
            BehResults.FV_Q2=bandpower(FVQuart2,sampFreq,[0.03,0.07]);
            BehResults.FV_Q3=bandpower(FVQuart3,sampFreq,[0.03,0.07]);
            BehResults.FV_Q4=bandpower(FVQuart4,sampFreq,[0.03,0.07]);
            BehResults.FVcontrol_Q1=bandpower(FVQuart1,sampFreq,[0.1,0.2]);
            BehResults.FVcontrol_Q2=bandpower(FVQuart2,sampFreq,[0.1,0.2]);
            BehResults.FVcontrol_Q3=bandpower(FVQuart3,sampFreq,[0.1,0.2]);
            BehResults.FVcontrol_Q4=bandpower(FVQuart4,sampFreq,[0.1,0.2]);            
            BehResults.FV=bandpower(FourthTC,sampFreq,[0.03,0.07]);
            BehResults.FVcontrol=bandpower(FourthTC,sampFreq,[0.1,0.2]);
            
            BehResults.Global_Q1=std(GlobalQuart1,0,1);
            BehResults.Global_Q2=std(GlobalQuart2,0,1);
            BehResults.Global_Q3=std(GlobalQuart3,0,1);
            BehResults.Global_Q4=std(GlobalQuart4,0,1);          
            BehResults.Global=std(globalSignal,0,1);
            BehResults.wm_Q1=std(wmQuart1,0,1);
            BehResults.wm_Q2=std(wmQuart2,0,1);
            BehResults.wm_Q3=std(wmQuart3,0,1);
            BehResults.wm_Q4=std(wmQuart4,0,1);          
            BehResults.wm=std(wmSignal,0,1);            
           
            for windNum = 1:size(ResultMatCSF,1)
                ResultMatCSF(windNum,[windNum:windNum+WindowSize-1]) = bandpower(FourthTC([windNum:windNum+WindowSize-1],1),sampFreq,[0.03,0.07]);
                ResultMatCSFcontrol(windNum,[windNum:windNum+WindowSize-1]) = bandpower(FourthTC([windNum:windNum+WindowSize-1],1),sampFreq,[0.1,0.2]);
                ResultMatGlobal(windNum,[windNum:windNum+WindowSize-1]) = std(globalSignal([windNum:windNum+WindowSize-1],1),0,1);
            end      
           % ResultMatCSF=zscore(nanmean(ResultMatCSF,1))';
           % ResultMatGlobal=zscore(nanmean(ResultMatGlobal,1))';
            ResultMatCSF=nanmean(ResultMatCSF,1)';
            ResultMatCSFcontrol=nanmean(ResultMatCSFcontrol,1)';
            ResultMatGlobal=nanmean(ResultMatGlobal,1)'; 
            BehResults.GlobalFVCorr=corr(ResultMatCSF,ResultMatGlobal);
            disp(corr(ResultMatCSF,ResultMatGlobal));
            beh_events.VigTC_global(1:endInd,1)=imresize(ResultMatGlobal,[endInd,1]);
            beh_events.VigTC_FV_Raw60vox(1:endInd,1)=imresize(ResultMatCSF,[endInd,1]);
            beh_events.VigTC_FVctr_Raw60vox(1:endInd,1)=imresize(ResultMatCSFcontrol,[endInd,1]);
            add_vars=struct2table(BehResults);
            beh_vars=[beh_vars,add_vars];
            save(LoadPath_behVars,'beh_vars');
            save(LoadPath_Events,'beh_events');      
        catch
            disp(['Error: ',ParcelNames{i,1}])
            continue
        end
    end        
    toc
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