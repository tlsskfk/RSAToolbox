function [AnalysisParameters] = ComputeParcellationReferenceRCA(fmriprep_table,ExperimentsDir,varargin)
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
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
[RSMsName] = VariableSetter('RSMsName',[],varargin);
%Subject or run level analysis. Will prompt request.
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);
[ParcelType] = VariableSetter('ParcelType',[],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','ReferenceRCA','TitleTextName','Select Analysis Parameters:');
            tempParcelNames=filePaths_RSMs.Properties.VariableNames;
        catch
            filePaths_RSMs=[];
            disp('No Existing Parameters!')
        end
        if ~isempty(filePaths_RSMs)
            for j = 1:height(filePaths_RSMs)
                try
                    load(filePaths_RSMs.(tempParcelNames{1,1}){j,1},'AnalysisParameters');
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
    SingleSelect=0; %Allows only a single value to be selected.
    [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect);    
    if ~isempty(EditParams)
        for i = 1:length(EditParams)
            AnalysisParameters.(EditParams{i,1})=[];
        end
    end
    SubjectOrRun=AnalysisParameters.SubjectOrRun;
    ParcelNames=AnalysisParameters.ParcelNames;   
    ParcelType=AnalysisParameters.ParcelType;  
else
    AnalysisParameters=struct;
end

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
    if length(unique(fmriprep_table.session))>1
        useIndicies=find(fmriprep_table.numRuns_bySes)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;        
    else
        useIndicies=find(fmriprep_table.numRuns_bySub)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    end
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end
AnalysisParameters.SubjectOrRun=SubjectOrRun;
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ReferenceRSMs';
AnalysisParameters.AnalysisType=AnalysisType;

if LoadActivationPatterns==1
    [filePaths_ActivationPattern,~,ActivationPatternName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ActivationPatterns','AnalysisName',ActivationPatternName,'TitleTextName','Select activation patterns for RSA:');
    filePaths_ActivationPattern=filePaths_ActivationPattern.WholeBrain;
end
AnalysisParameters.ActivationPatternName=ActivationPatternName;
%Allows you to set name for this particular analysis
if LoadActivationPatterns==1
    tempName=[strrep(ActivationPatternName,'ActivationPattern_','RSMs_'),gmAffix,ConfoundAffix,'_',SimType,UseStat];
else
    tempName=[AnalysisType,gmAffix,ConfoundAffix,'_',SimType,UseStat,afniSuffix];
end
if isempty(AnalysisName)
    tempName=strrep(tempName,'_cl_','_');
    tempName=strrep(tempName,'__','_');
    AnalysisName=uiEnterName(tempName,['Enter name for ',AnalysisType,newline,'analysis below:']);
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
    ParcelType=uiNameSelect([{'ParcelCoords','GroupParcel'}],'Select parcellations to run.');
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
else
    [ParcelNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','CoordParcels','GetAnalysisNames',1);
    ParcelNames=uiNameSelect([ParcelNames],'Select parcellations to run.');
    AllParcelTable=[];
    for i = 1:length(ParcelNames)
        [tempParcelFileNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','CoordParcels','AnalysisName',ParcelNames{i,1});
        AllParcelTable=[AllParcelTable,tempParcelFileNames];
    end
    ParcelNames=AllParcelTable.Properties.VariableNames(:);
end
   
AnalysisParameters.ParcelNames=ParcelNames;
numParcels=length(ParcelNames);
LoadSLName=0;
RunSearchlight=0;
if contains(ParcelNames,'Searchlight')
    RunSearchlight=1;
    if isempty(SearchlightName) && isempty(SearchlightParams)
        SingleSelect=1; %Allows only a single value to be selected.
        [LoadSLName] = uiNameSelect({'Yes','No'},'Load searchlight info?:',SingleSelect);
        if strcmpi(LoadSLName,'Yes')
            LoadSLName=1;
            SearchlightParams=1;
            [filePaths_SearchlightInfo,~,SearchlightName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SearchlightInfo','AnalysisName',SearchlightName,'TitleTextName','Select searchlight info:');
            filePaths_SearchlightInfo=filePaths_SearchlightInfo.(SearchlightName);
        else
            filePaths_SearchlightInfo=[];
            LoadSLName=0;         
            SingleSelect=1; %Allows only a single value to be selected.
            [slShape] = uiNameSelect({'Sphere','Cube'},'Select searchlight shape:',SingleSelect);
            slRadiusName=uiEnterName('3','Enter searchlight radius');
            slRadius=str2num(slRadiusName);   
            if strcmpi(slShape,'Sphere')
                maxVoxels=MaxSphereSLSize(slRadius,2);
            else
                maxVoxels=(slRadius*2+1)^3;
            end
            slVoxelThreshold=uiEnterName(num2str(round(maxVoxels*0.75)),['Enter searchlight radius ',newline,'Max # voxels: ',num2str(maxVoxels)]);
            SearchlightName=uiEnterName(['SearchlightInfo_Shape-',slShape,'_Rad-',slRadiusName,'_Thresh-',slVoxelThreshold],['Enter searchlight name:']);
            slVoxelThreshold=str2num(slVoxelThreshold);
            SearchlightParams.slShape=slShape;
            SearchlightParams.slRadius=slRadius;
            SearchlightParams.slVoxelThreshold=slVoxelThreshold;
        end
        
    end
else
    SearchlightName='';
end
SearchlightName=strrep(SearchlightName,'SearchlightInfo_','Searchlight_');
AnalysisParameters.SearchlightName=SearchlightName;
AnalysisParameters.SearchlightParams=SearchlightParams;
Parcels=cell(numParcels,1);
if strcmpi(ParcelType,'GroupParcel')
    for i = 1:numParcels
        if ~strcmpi(ParcelNames{i,1},'Searchlight')
            Parcels{i,1}=load(['Parcellations/',ParcelNames{i,1}],'UseMask','UseLabels');
        else
            Parcels{i,1}='Searchlight';
            ParcelNames{i,1}=SearchlightName;
        end
    end

    %Create mask that excludes voxels that aren't included in any parcellation.
    %Helpful to save time and minimize computation.
    analysisMask=[];

    if ~any(ismember(ParcelNames,SearchlightName)) && Compute_WM_RSM==0 && Compute_CSF_RSM==0 && Compute_GM_RSM==0
        for i = 1:numParcels
            UseMask=Parcels{i,1}.UseMask;
            if i ==1            
                analysisMask=single(UseMask>0);
            else
                try
                    analysisMask=analysisMask+single(UseMask>0);
                catch
                    disp('oops')
                    continue
                end
            end
            analysisMask=single(analysisMask>0);
        end
    end
else
   analysisMask=[];
end
BaseAnalysisMask=analysisMask;
if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end

AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.SearchlightParams=SearchlightParams;
AnalysisParameters.BaseTCName=BaseTCName;

%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
for dataInd=useIndicies
    %try
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
        descript1='desc-rsms'; %set file description name
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
        RunDur=cell(1);
        Use_ConfoundTCsByRun=cell(1);
        Use_ConfoundTCsBySubj=cell(1);
        Use_ActivationPatterns=[];
        Use_Events=cell(1);
        Use_TimeCourses=cell(1);
        Use_AM2Events=cell(1);
        Use_SearchlightInfo=[];
        Use_gmMask=[];    
        TR=fmriprep_table.TR(dataInd,1); 
        TrialNum=cell(1);     
        TrialOnsetTimes=cell(1);    
        %% load input data 
        % If by Subject, iterate through runs and place data in cell
        if LoadActivationPatterns==1
            LoadPath_ActivationPatterns=filePaths_ActivationPattern{dataInd,1};
            if strcmpi(UseStat,'T')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'tVals','brain_mask','DesignMatrix','RegressorNames','AfniInfo','AnalysisParameters');
                    Use_ActivationPatterns.ActVals=Use_ActivationPatterns.tVals;
                    Use_ActivationPatterns.tVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end
            elseif strcmpi(UseStat,'B')
                try
                    Use_ActivationPatterns = load(LoadPath_ActivationPatterns,'bVals','brain_mask','DesignMatrix','RegressorNames','AfniInfo','AnalysisParameters');
                    Use_ActivationPatterns.ActVals=Use_ActivationPatterns.bVals;
                    Use_ActivationPatterns.bVals=[];
                    brain_mask{1,1}=Use_ActivationPatterns.brain_mask;
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ActivationPatterns]);
                    continue 
                end 
            end
            brain_mask=repmat(brain_mask,[numRuns,1]);
        end   
        if LoadSLName==1
            LoadPath_SearchlightInfo=filePaths_SearchlightInfo{dataInd,1};
            Use_SearchlightInfo = load(LoadPath_SearchlightInfo,'SLcoords','SLinds','SLnumVoxels','AnalysisParameters');
            Use_SearchlightInfo.slShape=Use_SearchlightInfo.AnalysisParameters.slShape;
            Use_SearchlightInfo.slRadius=Use_SearchlightInfo.AnalysisParameters.slRadius;
            Use_SearchlightInfo.slVoxelThreshold=Use_SearchlightInfo.AnalysisParameters.slVoxelThreshold;
        else
            if RunSearchlight == 1
                Use_SearchlightInfo=SearchlightParams;
                Use_SearchlightInfo.SLcoords=[];
                Use_SearchlightInfo.SLinds=[];
                Use_SearchlightInfo.SLnumVoxels=[];
            end
        end
        for run=1:numRuns
            loadInd=dataInd+run-1;
            %skip previous errors
            if ~isempty(fmriprep_table.Error{loadInd,1})
                continue
            end    
            %% set load paths and variable names
            %pull load paths
            if LoadActivationPatterns==0 
                LoadPath_Events=filePaths_beh_events{loadInd,1};
                LoadPath_ConfoundTCs=filePaths_ConfoundTCs{loadInd,1};
                LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
                %Pull base timecourse data. If it doesn't exist, skip run.
                try
                    TempLoadData = load(LoadPath_BaseTC,'boldTCs','brain_mask');
                    BaseTC{count,1}=TempLoadData.boldTCs; 
                    brain_mask{count,1}=TempLoadData.brain_mask;
                    RunDur{count,1}=size(BaseTC{count,1},1);
                catch
                    disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
                    continue
                end 

                % pull confound TCs
                try
                    TempLoadData = load(LoadPath_ConfoundTCs,'ConfoundTCs');
                    if ~isempty(ConfoundNamesByRun)
                        ConfoundTCs=TempLoadData.ConfoundTCs(:,TempLoadData.ConfoundTCs.Properties.VariableNames(ismember(TempLoadData.ConfoundTCs.Properties.VariableNames,ConfoundNamesByRun)));
                        Use_ConfoundTCsByRun{count,1}=table2array(ConfoundTCs);
                        Use_ConfoundTCsByRun{count,2}=ConfoundTCs.Properties.VariableNames(:);
                    end
                    if ~isempty(ConfoundNamesBySubj)
                        ConfoundTCs=TempLoadData.ConfoundTCs(:,TempLoadData.ConfoundTCs.Properties.VariableNames(ismember(TempLoadData.ConfoundTCs.Properties.VariableNames,ConfoundNamesBySubj)));
                        Use_ConfoundTCsBySubj{count,1}=table2array(ConfoundTCs);
                        Use_ConfoundTCsBySubj{count,2}=ConfoundTCs.Properties.VariableNames(:);
                    end            
                catch
                    disp(['No confound regressors-- input file or variable doesnt exist: ',LoadPath_ConfoundTCs]);
                end   

                %Pull experiment- or  behavior-based events, timecourses and AM2 events
                if ~isempty(LoadPath_Events)
                    try
                        TempLoadData = load(LoadPath_Events,'beh_events');
                        if ~isempty(EventNames)
                            Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,EventNames)));
                            Use_Events{count,1}=table2array(Events);
                            Use_Events{count,2}=Events.Properties.VariableNames(:);
                        end 
                        if ~isempty(TimeCourseNames)
                            TimeCourses=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,TimeCourseNames)));
                            Use_TimeCourses{count,1}=table2array(TimeCourses);
                            Use_TimeCourses{count,2}=TimeCourses.Properties.VariableNames(:);
                        end 
                        if ~isempty(AM2Event_Names)
                            AM2Events=TempLoadData.beh_events(:,TempLoadData.beh_events.Properties.VariableNames(ismember(TempLoadData.beh_events.Properties.VariableNames,AM2Event_Names)));
                            Use_AM2Events{count,1}=table2array(AM2Events);
                            Use_AM2Events{count,2}=AM2Events.Properties.VariableNames(:);
                        end                
                    catch
                        disp(['No events-- input file or variable doesnt exist: ',LoadPath_Events]);
                    end 
                    try
                        TempLoadData = load(LoadPath_Events,'beh_events');     
                        TrialNum{count,1}=TempLoadData.beh_events.TrialNum; 
                        TrialOnsetTimes{count,1}=TempLoadData.beh_events.TrialOnsetTime; 
                    catch
                       disp(['Event load error-- skipping',LoadPath_Events]); 
                       continue
                    end                  
                end
            end
            if ~isempty(gmMask) && isempty(Use_gmMask)
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
            end

            if Compute_WM_RSM == 1
                try
                    TempLoadData = load(filePaths_gmMask{loadInd,1},'WM_probseg');
                    Use_wmMask=single(TempLoadData.WM_probseg > 0.99);
                catch
                    disp(['White matter mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                    Use_wmMask=[];
                end
            else
                Use_wmMask=[];
            end 

            if Compute_CSF_RSM == 1
                try
                    TempLoadData = load(filePaths_gmMask{loadInd,1},'CSF_probseg');
                    Use_csfMask=single(TempLoadData.CSF_probseg > 0.9);
                catch
                    disp(['CSF mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
                    Use_csfMask=[];
                end
            else
                Use_csfMask=[];
            end     
            if Compute_GM_RSM == 1
                try
                    TempLoadData = load(filePaths_gmMask{loadInd,1},'GM_probseg');
                    Use_invgmMask=single(TempLoadData.GM_probseg < 0.01);
                catch
                    continue
                    disp(['invGM mask error-- mask not applied']);
                    continue
                    Use_invgmMask=[];
                end
            else
                Use_invgmMask=[];
            end              
            count=count+1;
        end  

        if count==1
            disp(['Skipping subject-- no input files or variables exist']);
            continue
        end  
        if Compute_CSF_RSM == 0 && Compute_WM_RSM == 0 && Compute_GM_RSM == 0
            if isempty(analysisMask)
                analysisMask=Use_gmMask;
            else
                analysisMask=analysisMask.*Use_gmMask;
            end
        else
            analysisMask=(brain_mask{1,1}*0)+1;        
        end
        if strcmpi(ParcelType,'ParcelCoords')
            tempParcelPaths=table2cell(AllParcelTable(dataInd,:));
            Parcels=cell(size(tempParcelPaths,2),1);
            for j = 1:size(tempParcelPaths,2)
                Parcels{j,1}=load(tempParcelPaths{1,j},'UseMask','UseLabels');          
                if unique(Parcels{j,1}.UseMask(:)) == 0
                    disp(['Skipping ',fmriprep_table.sub{dataInd,1},'-- No parcel data']);
                    continue
                end
            end
        end
        
        %% Run analysis here!!
        [All_RSMs,All_rsm_masks,All_ConfoundRSMs,AnalysisParameters.RegressorNames,AnalysisParameters.AfniInfo,SearchlightResults] = GetParcellationRSMs(...
            Use_Events,Use_TimeCourses,Use_ConfoundTCsByRun,...
            BaseTC,brain_mask,ConditionNames,Parcels,...
            'SearchlightInfo',Use_SearchlightInfo,...
            'ActivationPatterns',Use_ActivationPatterns,...
            'UseAfni',UseAfni,...
            'Use3DReml',Use3DReml,...
            'TrialOnsetTimes',TrialOnsetTimes,...
            'AfniWorkDir',AfniWorkDir,...
            'TR',TR,...
            'TrialNum',TrialNum,...            
            'csfMask',Use_csfMask,...
            'wmMask',Use_wmMask,...
            'invgmMask',Use_invgmMask,...
            'Compute_MeanAct_RSM',Compute_MeanAct_RSM,...
            'gmMask',Use_gmMask,...
            'SimType',SimType,...
            'OutputConfoundRSMs',0,...
            'parGLM',1,...
            'Compute_CondReg_CorrMat',Compute_CondReg_CorrMat,...
            'ContrastNames',cNames,...
            'Normalize','zscore',...
            'ConfoundTCsBySubj',Use_ConfoundTCsBySubj,...
            'ResampleSizes',RunDur,...
            'AM2',AM2,...
            'AM2_Events',Use_AM2Events,...
            'TimecourseShift',0,...
            'ParcellationMask',analysisMask,...
            'ResampleSizesBoldTC',[],...
            'NormWithinRun',[],...
            'NormAcrossRun','zscore',...
            'Contrasts',cMat,...
            'BatchSize',BatchSize,...
            'ResampleToFit','Y',...
            'UseStat',UseStat);   
        for i = 1:numParcels
            if strcmpi(ParcelNames{i,1},'Searchlight')
                RSMs=All_RSMs{i,1};
                rsm_mask=All_rsm_masks{i,1};                
                save(SaveNames{i,1},'RSMs','rsm_mask','SearchlightResults','AnalysisParameters');
            else 
                RSMs=All_RSMs{i,1};
                rsm_mask=All_rsm_masks{i,1};
                if ~isempty(All_ConfoundRSMs)
                    confoundRSMs=All_ConfoundRSMs{i,1};
                    save(SaveNames{i,1},'RSMs','rsm_mask','confoundRSMs','AnalysisParameters');   
                else
                    save(SaveNames{i,1},'RSMs','rsm_mask','AnalysisParameters');   
                end
            end
        end 
        toc
%     catch
%        disp('Error!');
%     end
end 
end