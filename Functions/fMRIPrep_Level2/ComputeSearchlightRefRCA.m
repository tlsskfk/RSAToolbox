function [AnalysisParameters] = ComputeSearchlightRefRCA(fmriprep_table,ExperimentsDir,varargin)
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
[ParcelNames] = VariableSetter('ParcelNames',[],varargin);
[RSMsName] = VariableSetter('RSMsName',[],varargin);
[RefInfo] = VariableSetter('RefInfo',[],varargin);
%Subject or run level analysis. Will prompt request.
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);
CleanRSMs = VariableSetter('CleanRSMs',['No'],varargin);
BatchSize = VariableSetter('BatchSize',1000,varargin);
nThreshold = VariableSetter('nThreshold',[],varargin);
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
    RSMsName=AnalysisParameters.RSMsName;    
    ParcelNames=AnalysisParameters.ParcelNames;    
    RefInfo=AnalysisParameters.RefInfo;  
    FigureInfo=AnalysisParameters.FigureInfo;
    CleanRSMs=AnalysisParameters.CleanRSMs;
else
    AnalysisParameters=struct;
end

if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end

ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');

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
AnalysisParameters.SubjectOrRun=SubjectOrRun;

if isempty(RefInfo)
    [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','GroupType','GroupAnalysis','FileNames','slRSMs.mat');
else
    [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','FileNames','slRSMs.mat','GroupType','GroupAnalysis',...
        'AnalysisNames',RefInfo.AnalysisNames,'TableNames',RefInfo.TableNames,'ParcellationNames',RefInfo.ParcellationNames,'LoadVarNames_Parcel',RefInfo.LoadVarNames_Parcel);
end
RefParcelNames=RefInfo.ParcellationNames;
RefAppend=[RefInfo.TableNames{1,1},'_',RefInfo.AnalysisNames{1,1}];
RefAppend=strrep(RefAppend,'_RSMs','');
load(AllPaths_RefRSMs{1,1});
if isempty(nThreshold)
     nThreshold=uiEnterName('',['Enter refRSM overlap theshold out of ', num2str(max(OverlapVec))]);
     nThreshold=str2num(nThreshold);
end
OverlapMap=single(slMask)*0;
OverlapMap(slMask)=OverlapVec;
OverlapMap=OverlapMap >= nThreshold;
RefVec=OverlapMap(:);
SelectInd = OverlapVec >= nThreshold;
slRSMs=slRSMs(:,SelectInd);
slZRSMs=slZRSMs(:,SelectInd);
if isempty(FigureInfo)
    SingleSelect=0; %Allows only a single value to be selected.
    [FigureInfo] = uiNameSelect({'None','MakeFigures','ComputeReliability','MakeBrainMaps'},'Select summary figures to make: ',SingleSelect);    
end
MakeFigs=0;
ComputeReliability=0;
MakeBrainMaps=0;
if any(contains(FigureInfo,'MakeFigures'))
    MakeFigs=1;
end
if any(contains(FigureInfo,'ComputeReliability'))
    ComputeReliability=1;
end
if any(contains(FigureInfo,'MakeBrainMaps'))
    MakeBrainMaps=1;
end
%% Compile filepaths for input files for the analysis
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',RSMsName,'TitleTextName',['Select RSMs for',newline,'Representational Connectivity Analysis:']);
if isempty(ParcelNames)
    ParcelNames=filePaths_RSMs.Properties.VariableNames(:);
    %ParcelNames=ParcelNames(ismember(ParcelNames,RefParcelNames),:);
    ParcelNames=uiNameSelect(ParcelNames,'Select Searchlight to run.',1);
else
    ParcelNames=ParcelNames(ismember(ParcelNames,filePaths_RSMs.Properties.VariableNames(:)),:);
    ParcelNames=ParcelNames(ismember(ParcelNames,RefParcelNames),:);
end
filePaths_RSMs=filePaths_RSMs(:,ParcelNames);
numParcels=length(ParcelNames);
if isempty(CleanRSMs)
    CleanRSMs=uiNameSelect({'Yes','No'},'Use cleaned RSMs?');
end
if strcmpi(CleanRSMs,'Yes')
    clAppend='_cl';
else
    clAppend='';
end

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='RefRCA';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName([strrep(RSMsName,'RSMs_',[AnalysisType,'_']),'_',RefAppend,clAppend],['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

if Overwrite==1 || ~isfield(AnalysisParameters,'RunDate')
    AnalysisParameters.RunDate=genDateString;
end
AnalysisParameters.AnalysisType=AnalysisType;
AnalysisParameters.AnalysisName=AnalysisName;
AnalysisParameters.Parcellation=ParcelNames;
AnalysisParameters.RSMsName=RSMsName;
AnalysisParameters.RefInfo=RefInfo;
AnalysisParameters.FigureInfo=FigureInfo;
AnalysisParameters.CleanRSMs=CleanRSMs;
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
        runNum=fmriprep_table.run(dataInd,1);
        SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
        SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
    end
    descript1='desc-RefRCA'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    
    for parcelNum=1:numParcels
        SaveNames=cell(1); 
        try
            ParcelName=ParcelNames{parcelNum,1};
        catch
            ParcelName=ParcelNames;
        end
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
            LoadPath_RSMs=filePaths_RSMs.(ParcelName){dataInd,1};
        catch
            LoadPath_RSMs=[];
        end
        %Pull base timecourse data. If it doesn't exist, skip run.
        if isempty(LoadPath_RSMs)
            continue
        end
        try
            TempLoadData = load(LoadPath_RSMs,'RSMs','rsm_mask');
            RSMs=single(TempLoadData.RSMs);
            rsm_mask=single(TempLoadData.rsm_mask);
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_RSMs]);
            continue
        end               
        t=rsm_mask(:);
        rsm_mask=(single(rsm_mask>0)+single(OverlapMap>0))==2;
        RSMs=RSMs(:,RefVec(t>0)==1);
        tempRefRSMs=slRSMs(:,t(RefVec>0)==1);
        tempRefZRSMs=slZRSMs(:,t(RefVec>0)==1);

        [Ref_RF] = BatchRefRFCorr(RSMs,tempRefRSMs,BatchSize);
%         RCAMat=corr(RSMs,tempRefRSMs);
%         Ref_RF=single(RCAMat(eye(length(RCAMat))==1)');
%         Ref_RC=single((RCAMat+permute(RCAMat,[2,1]))/2);
%         Ref_RC(eye(length(RCAMat))==1)=nan;    
% 
%         RCAMat=corr(RSMs,tempRefZRSMs);
%         RefZ_RF=single(RCAMat(eye(length(RCAMat))==1)');
        [RefZ_RF] = BatchRefRFCorr(RSMs,tempRefZRSMs,BatchSize);
        RFMap=single(rsm_mask*0);
        RFzMap=RFMap;
        RFMap(rsm_mask==1)=Ref_RF;
        RFzMap(rsm_mask==1)=RefZ_RF;
%         RefZ_RC=single((RCAMat+permute(RCAMat,[2,1]))/2);
%         RefZ_RC(eye(length(RCAMat))==1)=nan; 
        save(SaveNames{1,1},'RFMap','RFzMap');
        %save(SaveNames{1,1},'rsm_mask','Ref_RF','Ref_RC','RefZ_RF','RefZ_RC');
    end     
toc    
end
try
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames{1,1},'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames{1,1},'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames{1,1},'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames{1,1},'/'];
    end
catch
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames,'/'];
    end    
end    
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);
    mkdir(GroupBrainMapsDir);
end
try
    [RFMaps] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','RFMap',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','RefRCA',...
        'AnalysisName',AnalysisName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelNames{1,1});
catch
    [RFMaps] = CompileND(ExperimentsDir,fmriprep_table,...
        'LoadVarName','RFMap',...
        'LoadVarFormat','Array',...
        'TableOrArray','Array',...
        'DataType','ByParcellation',...
        'AnalysisType','RefRCA',...
        'AnalysisName',AnalysisName,...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelNames);    
end    
ThreshMask=single(sum(single(~isnan(RFMaps)),4)>round(0.75*size(RFMaps,4)));
ThreshMask=repmat(ThreshMask,[1,1,1,size(RFMaps,4)]);
RFMaps(ThreshMask==0)=nan;
[groupRefRF_t,groupRefRF_z,groupRefRF_mean]=getTval(atanh(RFMaps),4);
groupRefRF_t(isinf(groupRefRF_t))=nan;
groupRefRF_z(isinf(groupRefRF_z))=nan;
groupRefRF_mean(isinf(groupRefRF_mean))=nan;
BrainMap=cat(4,groupRefRF_mean,groupRefRF_t,groupRefRF_z);
SaveBrik_3mmMNI(BrainMap,{'group_mean','group_t','group_z'},[GroupBrainMapsDir,'groupRefRF_SearchlightMaps']);
save([GroupAnalysisDir,'groupRefRF_SearchlightMaps'],'groupRefRF_mean','groupRefRF_t','groupRefRF_z','AnalysisParameters');    
end

