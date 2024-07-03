function [AnalysisParameters] = ComputeSearchlight2SplitVarRefRCA(fmriprep_table,ExperimentsDir,varargin)
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
SplitRef = VariableSetter('SplitRef',[],varargin);
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
    SplitRef=AnalysisParameters.SplitRef;
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

if isempty(SplitRef)
    SingleSelect=1; %Allows only a single value to be selected.
    [SplitRef] = uiNameSelect({'Yes','No'},'Use split var reference RSMs:',SingleSelect);
end

if strcmpi(SplitRef,'Yes')
    UseSplitRef=1;
    if isempty(RefInfo)
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','SplitVarRSMs','GroupType','GroupAnalysis','FileNames','slRSMs.mat');
    else
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','SplitVarRSMs','FileNames','slRSMs.mat','GroupType','GroupAnalysis',...
            'AnalysisNames',RefInfo.AnalysisNames,'TableNames',RefInfo.TableNames,'ParcellationNames',RefInfo.ParcellationNames,'LoadVarNames_Parcel',RefInfo.LoadVarNames_Parcel);
    end
    UseSplitRef=['_SplitRef'];
else
    UseSplitRef=0;
    if isempty(RefInfo)
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','GroupType','GroupAnalysis','FileNames','slRSMs.mat');
    else
        [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','FileNames','slRSMs.mat','GroupType','GroupAnalysis',...
            'AnalysisNames',RefInfo.AnalysisNames,'TableNames',RefInfo.TableNames,'ParcellationNames',RefInfo.ParcellationNames,'LoadVarNames_Parcel',RefInfo.LoadVarNames_Parcel);
    end
    SplitRefAppend=[];
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
if UseSplitRef==0
    Split1_slRSMs=slRSMs;
    Split2_slRSMs=slRSMs;
    Split1_slzRSMs=slZRSMs;
    Split2_slzRSMs=slZRSMs;    
end
Split1_slRSMs=Split1_slRSMs(:,SelectInd);
Split1_slZRSMs=Split1_slZRSMs(:,SelectInd);
Split2_slRSMs=Split2_slRSMs(:,SelectInd);
Split2_slZRSMs=Split2_slZRSMs(:,SelectInd);
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
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SplitVarRSMs','AnalysisName',RSMsName,'TitleTextName',['Select RSMs for',newline,'Representational Connectivity Analysis:']);
if isempty(ParcelNames)
    ParcelNames=filePaths_RSMs.Properties.VariableNames(:);
    ParcelNames=ParcelNames(ismember(ParcelNames,RefParcelNames),:);
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
AnalysisType='SplitVarRefRCA';
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
AnalysisParameters.SplitRef=SplitRef;
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
            TempLoadData = load(LoadPath_RSMs,'Split1_RSMs','Split2_RSMs','brain_mask');
            RSMs=single(TempLoadData.RSMs);
            rsm_mask=single(TempLoadData.rsm_mask);
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_RSMs]);
            continue
        end               
        t=brain_mask(:);
        brain_mask=(single(brain_mask>0)+single(OverlapMap>0))==2;
        
        Split1_RSMs=Split1_RSMs(:,RefVec(t>0)==1);
        tempRefRSMs=Split1_slRSMs(:,t(RefVec>0)==1);
        tempRefZRSMs=Split1_slZRSMs(:,t(RefVec>0)==1);
        
        RCAMat=corr(Split1_RSMs,tempRefRSMs);
        Ref_RF=single(RCAMat(eye(length(RCAMat))==1)');
%         Ref_RC=single((RCAMat+permute(RCAMat,[2,1]))/2);
%         Ref_RC(eye(length(RCAMat))==1)=nan;    

        RCAMat=corr(Split1_RSMs,tempRefZRSMs);
        RefZ_RF=single(RCAMat(eye(length(RCAMat))==1)');
        
        Split1_RFMap=single(brain_mask*0);
        Split1_RFzMap=Split1_RFMap;
        Split1_RFMap(brain_mask==1)=Ref_RF;
        Split1_RFzMap(brain_mask==1)=RefZ_RF;
        
        Split2_RSMs=Split2_RSMs(:,RefVec(t>0)==1);
        tempRefRSMs=Split2_slRSMs(:,t(RefVec>0)==1);
        tempRefZRSMs=Split2_slZRSMs(:,t(RefVec>0)==1);
        
        RCAMat=corr(RSMs,tempRefRSMs);
        Ref_RF=single(RCAMat(eye(length(RCAMat))==1)');
%         Ref_RC=single((RCAMat+permute(RCAMat,[2,1]))/2);
%         Ref_RC(eye(length(RCAMat))==1)=nan;    

        RCAMat=corr(RSMs,tempRefZRSMs);
        RefZ_RF=single(RCAMat(eye(length(RCAMat))==1)');
        
        Split2_RFMap=single(brain_mask*0);
        Split2_RFzMap=RFMap;
        Split2_RFMap(brain_mask==1)=Ref_RF;
        Split2_RFzMap(brain_mask==1)=RefZ_RF;        
%         RefZ_RC=single((RCAMat+permute(RCAMat,[2,1]))/2);
%         RefZ_RC(eye(length(RCAMat))==1)=nan; 
        SplitDiff_RFMap=Split2_RFMap-Split1_RFMap;
        SplitDiff_RFzMap=Split2_RFzMap-Split1_RFzMap;
        save(SaveNames{1,1},'Split1_RFMap','Split1_RFzMap','Split2_RFMap','Split2_RFzMap','SplitDiff_RFMap','SplitDiff_RFzMap');
        %save(SaveNames{1,1},'rsm_mask','Ref_RF','Ref_RC','RefZ_RF','RefZ_RC');
    end     
toc    
end
if bySS == 1
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames{1,1},'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelNames{1,1},'/'];
else
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames{1,1},'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelNames{1,1},'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);
    mkdir(GroupBrainMapsDir);
end
[RFMaps] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','Split1_RFMap',...
    'LoadVarFormat','Array',...
    'TableOrArray','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','SplitVarRefRCA',...
    'AnalysisName',AnalysisName,...
    'SubjectOrRun',SubjectOrRun,...
    'ParcelName',ParcelNames{1,1});
ThreshMask=single(sum(single(~isnan(RFMaps)),4)>round(0.75*size(RFMaps,4)));
ThreshMask=repmat(ThreshMask,[1,1,1,size(RFMaps,4)]);
RFMaps(ThreshMask==0)=nan;
[Split1_groupRefRF_t,Split1_groupRefRF_z,Split1_groupRefRF_mean]=getTval(atanh(RFMaps),4);
Split1_groupRefRF_t(isinf(Split1_groupRefRF_t))=nan;
Split1_groupRefRF_z(isinf(Split1_groupRefRF_z))=nan;
Split1_groupRefRF_mean(isinf(Split1_groupRefRF_mean))=nan;
BrainMap=cat(4,Split1_groupRefRF_mean,Split1_groupRefRF_t,Split1_groupRefRF_z);

[RFMaps] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','Split2_RFMap',...
    'LoadVarFormat','Array',...
    'TableOrArray','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','SplitVarRefRCA',...
    'AnalysisName',AnalysisName,...
    'SubjectOrRun',SubjectOrRun,...
    'ParcelName',ParcelNames{1,1});
ThreshMask=single(sum(single(~isnan(RFMaps)),4)>round(0.75*size(RFMaps,4)));
ThreshMask=repmat(ThreshMask,[1,1,1,size(RFMaps,4)]);
RFMaps(ThreshMask==0)=nan;
[Split2_groupRefRF_t,Split2_groupRefRF_z,Split2_groupRefRF_mean]=getTval(atanh(RFMaps),4);
Split2_groupRefRF_t(isinf(Split2_groupRefRF_t))=nan;
Split2_groupRefRF_z(isinf(Split2_groupRefRF_z))=nan;
Split2_groupRefRF_mean(isinf(Split2_groupRefRF_mean))=nan;
BrainMap=cat(4,BrainMap,Split2_groupRefRF_mean,Split2_,groupRefRF_t,Split2_groupRefRF_z);

[RFMaps] = CompileND(ExperimentsDir,fmriprep_table,...
    'LoadVarName','SplitDiff_RFMap',...
    'LoadVarFormat','Array',...
    'TableOrArray','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','SplitVarRefRCA',...
    'AnalysisName',AnalysisName,...
    'SubjectOrRun',SubjectOrRun,...
    'ParcelName',ParcelNames{1,1});
ThreshMask=single(sum(single(~isnan(RFMaps)),4)>round(0.75*size(RFMaps,4)));
ThreshMask=repmat(ThreshMask,[1,1,1,size(RFMaps,4)]);
RFMaps(ThreshMask==0)=nan;
[SplitDiff_groupRefRF_t,SplitDiff_groupRefRF_z,SplitDiff_groupRefRF_mean]=getTval(atanh(RFMaps),4);
SplitDiff_groupRefRF_t(isinf(SplitDiff_groupRefRF_t))=nan;
SplitDiff_groupRefRF_z(isinf(SplitDiff_groupRefRF_z))=nan;
SplitDiff_groupRefRF_mean(isinf(SplitDiff_groupRefRF_mean))=nan;
BrainMap=cat(4,BrainMap,SplitDiff_groupRefRF_mean,SplitDiff_groupRefRF_t,SplitDiff_groupRefRF_z);
SaveBrik_3mmMNI(BrainMap,{'Split1_group_mean','Split1_group_t','Split1_group_z','Split2_group_mean','Split2_group_t','Split2_group_z','SplitDiff_group_mean','SplitDiff_group_t','SplitDiff_group_z'},[GroupBrainMapsDir,'SplitVar_groupRefRF_SearchlightMaps']);

save([GroupAnalysisDir,'SplitVargroupRefRF_SearchlightMaps'],'Split1_groupRefRF_mean','Split1_groupRefRF_t','Split1_groupRefRF_z','Split2_groupRefRF_mean','Split2_groupRefRF_t','Split2_groupRefRF_z','SplitDiff_groupRefRF_mean','SplitDiff_groupRefRF_t','SplitDiff_groupRefRF_z','AnalysisParameters');    
end

