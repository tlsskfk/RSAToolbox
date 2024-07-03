function [AnalysisParameters] = ComputeParcellationRefRCA(fmriprep_table,ExperimentsDir,varargin)
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
CleanRSMs = VariableSetter('CleanRSMs',[],varargin);
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
tableVars=fmriprep_table.Properties.VariableNames;
SubOrRunOpts={'Subject','Run'};
if any(ismember(tableVars,'numRuns_bySes'))
    SubOrRunOpts=[SubOrRunOpts,{'Session'}];
end
if any(ismember(tableVars,'numRuns_byGroup'))
    SubOrRunOpts=[SubOrRunOpts,{'Group'}];
end
if ~any(ismember(tableVars,'numRuns_bySub'))
    fmriprep_table.numRuns_bySub=fmriprep_table.numRuns;
end
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect(SubOrRunOpts,'Perform analysis by subject or by run:',SingleSelect);
end
bySes=0;
bySS=0;
bySub=0;
byGroup=0;
byRun=0;
groupName=[];
ssAppend=[];
if strcmpi(SubjectOrRun,'Subject')
    bySub=1;
    useIndicies=find(fmriprep_table.numRuns_bySub)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    ssAppend='_bySub';
elseif strcmpi(SubjectOrRun,'Group')
    byGroup=1;
    useIndicies=find(fmriprep_table.numRuns_byGroup)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_byGroup; 
    groupName=uiEnterName('','Enter group name.');
    ssAppend=['_by',groupName];
elseif strcmpi(SubjectOrRun,'Session')
    bySes=1;
    useIndicies=find(fmriprep_table.numRuns_bySes)';
    fmriprep_table.numRuns=fmriprep_table.numRuns_bySes; 
    ssAppend='_bySes';
else    
    byRun=1;
    useIndicies=[1:TotalRuns];
    ssApend='_byRun';
end

if isempty(RefInfo)
    [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','GroupType','GroupAnalysis','FileNames','ParcelRSMs.mat');
else
    [AllPaths_RefRSMs,~,~,~,RefInfo] = GroupDirSearch(ExperimentsDir,'AnalysisTypes','RSMs','FileNames','ParcelRSMs.mat','GroupType','GroupAnalysis',...
        'AnalysisNames',RefInfo.AnalysisNames,'TableNames',RefInfo.TableNames,'ParcellationNames',RefInfo.ParcellationNames,'LoadVarNames_Parcel',RefInfo.LoadVarNames_Parcel);
end
RefParcelNames=RefInfo.ParcellationNames;
RefAppend=[RefInfo.TableNames{1,1},'_',RefInfo.AnalysisNames{1,1}];
RefAppend=strrep(RefAppend,'_RSMs','');
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
    ParcelNames=ParcelNames(ismember(ParcelNames,RefParcelNames),:);
    ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
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
    grpAppend=[];
    if bySes == 1
        if fmriprep_table.numRuns_bySes(dataInd,1)>1 
            runNum=fmriprep_table.run(dataInd,1);
            SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
            grpAppend='_grp-bySes';
        end    
    
    elseif bySub==1
        if fmriprep_table.numRuns_bySub(dataInd,1)>1
            runNum=fmriprep_table.run(dataInd,1);
            sesNum=fmriprep_table.session(dataInd,1);
            SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_ses-0',num2str(sesNum)],'');
            SaveName=strrep(SaveName,['_ses-',num2str(sesNum)],''); 
            grpAppend='_grp-bySub';
        end    
    elseif byGroup==1
        if fmriprep_table.numGroup(dataInd,1)>1
            runNum=fmriprep_table.run(dataInd,1);
            sesNum=fmriprep_table.session(dataInd,1);   
            groupNum=fmriprep_table.groupNum(dataInd,1);  
            SaveName=strrep(SaveName,['_run-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_run-',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_ses-0',num2str(runNum)],'');
            SaveName=strrep(SaveName,['_ses-',num2str(runNum)],''); 
            SaveName=[SaveName,ssAppend,num2str(groupNum)];
            grpAppend=['_grp-by',groupName];
        end
    end
    if isempty(grpAppend)
        grpAppend='_grp-byRun';
    end
    descript1=['desc-RefRCA',grpAppend]; %set file description name
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
        if bySub==1
            numRuns=fmriprep_table.numRuns_bySub(dataInd,1);
        elseif bySes==1
            numRuns=fmriprep_table.numRuns_bySes(dataInd,1);
        elseif byGroup==1
            numRuns=fmriprep_table.numRuns_byGroup(dataInd,1);
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
            TempLoadData = load(LoadPath_RSMs,'RSMs');
            RSMs=single(TempLoadData.RSMs); 
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_RSMs]);
            continue
        end               
        %% Run analysis here!!
        parcelName=ParcelName;
        
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
        if strcmpi(CleanRSMs,'Yes')
            AllRSMs=load(AllPaths_RefRSMs{find(ismember(RefParcelNames,ParcelName)),1},'ParcelRSMs_p50','ParcelZRSMs_p50');
            AllRSMs.ParcelRSMs=AllRSMs.ParcelRSMs_p50;
            AllRSMs.ParcelZRSMs=AllRSMs.ParcelZRSMs_p50;
        else
            AllRSMs=load(AllPaths_RefRSMs{find(ismember(RefParcelNames,ParcelName)),1},'ParcelRSMs','ParcelZRSMs');
        end

        if size(AllRSMs.ParcelRSMs,2)==1
            Ref_RF=array2table(corr(table2array(AllRSMs.ParcelRSMs),RSMs),'VariableNames',UseLabels);
            Ref_RF_Z = Ref_RF;
            RefZ_RF=array2table(corr(table2array(AllRSMs.ParcelZRSMs),RSMs),'VariableNames',UseLabels);
            RefZ_RF_Z = RefZ_RF;  
            Ref_RC=[];
            RefZ_RC=[];
            Ref_RC_WholeCM=[];
            Ref_RC_Degree=[];
            Ref_RC_Z=[];
            RefZ_RC_WholeCM=[];
            RefZ_RC_Degree=[];
            RefZ_RC_Z=[];
            save(SaveNames{1,1},'Ref_RF','Ref_RC','RefZ_RF','RefZ_RC','Ref_RC_WholeCM','Ref_RC_Degree','Ref_RC_Z','Ref_RF_Z','RefZ_RC_WholeCM','RefZ_RC_Degree','RefZ_RC_Z','RefZ_RF_Z');
        else    
            RCAMat=corr(table2array(AllRSMs.ParcelRSMs),RSMs);
            Ref_RF=array2table(RCAMat(eye(length(RCAMat))==1)','VariableNames',UseLabels);
            FullMat=(RCAMat+permute(RCAMat,[2,1]))/2;
            FullMat(eye(length(RCAMat))==1)=nan;
            Ref_RC=array2table(FullMat,'RowNames',UseLabels,'VariableNames',UseLabels);

            [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(Ref_RC),'ThresholdType','rank','Threshold',0.5);
            rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
            for nameNum=1:length(rowNames50)
                rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
            end
            tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);

            [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(Ref_RC),'ThresholdType','rank','Threshold',0.25);
            rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
            for nameNum=1:length(rowNames25)
                rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
            end        
            tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);

            [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(Ref_RC),'ThresholdType','rank','Threshold',0.10);
            rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
            for nameNum=1:length(rowNames10)
                rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
            end        
            tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);

            [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(Ref_RC),'ThresholdType','rank','Threshold',0.05);
            rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
            for nameNum=1:length(rowNames05)
                rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
            end        
            tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);

            Ref_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
            Ref_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
            Ref_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
            Ref_RC_Z=array2table(vertRSM2SymRSM(nan_zscore(mat2uppertriuvectormat(table2array(Ref_RC))),[],1),'VariableNames',Ref_RC.Properties.VariableNames,'RowNames',Ref_RC.Properties.RowNames);
            Ref_RF_Z = array2table(nan_zscore(table2array(Ref_RF)')','VariableNames',Ref_RF.Properties.VariableNames);        


            RCAMat=corr(table2array(AllRSMs.ParcelZRSMs),RSMs);
            RefZ_RF=array2table(RCAMat(eye(length(RCAMat))==1)','VariableNames',UseLabels);
            FullMat=(RCAMat+permute(RCAMat,[2,1]))/2;
            FullMat(eye(length(RCAMat))==1)=nan;
            RefZ_RC=array2table(FullMat,'RowNames',UseLabels,'VariableNames',UseLabels);

            [tempResultsByNodeRank50,tempResultsByCMRank50] = ConnectivityMatrixStats(table2array(RefZ_RC),'ThresholdType','rank','Threshold',0.5);
            rowNames50=tempResultsByNodeRank50.Properties.VariableNames';
            for nameNum=1:length(rowNames50)
                rowNames50{nameNum,1}=[rowNames50{nameNum,1},'_rank50'];
            end
            tempResultsByNodeRank50=array2table(table2array(tempResultsByNodeRank50)','VariableNames',UseLabels,'RowNames',rowNames50);

            [tempResultsByNodeRank25,tempResultsByCMRank25] = ConnectivityMatrixStats(table2array(RefZ_RC),'ThresholdType','rank','Threshold',0.25);
            rowNames25=tempResultsByNodeRank25.Properties.VariableNames';
            for nameNum=1:length(rowNames25)
                rowNames25{nameNum,1}=[rowNames25{nameNum,1},'_rank25'];
            end        
            tempResultsByNodeRank25=array2table(table2array(tempResultsByNodeRank25)','VariableNames',UseLabels,'RowNames',rowNames25);

            [tempResultsByNodeRank10,tempResultsByCMRank10] = ConnectivityMatrixStats(table2array(RefZ_RC),'ThresholdType','rank','Threshold',0.10);
            rowNames10=tempResultsByNodeRank10.Properties.VariableNames';
            for nameNum=1:length(rowNames10)
                rowNames10{nameNum,1}=[rowNames10{nameNum,1},'_rank10'];
            end        
            tempResultsByNodeRank10=array2table(table2array(tempResultsByNodeRank10)','VariableNames',UseLabels,'RowNames',rowNames10);

            [tempResultsByNodeRank05,tempResultsByCMRank05] = ConnectivityMatrixStats(table2array(RefZ_RC),'ThresholdType','rank','Threshold',0.05);
            rowNames05=tempResultsByNodeRank05.Properties.VariableNames';
            for nameNum=1:length(rowNames05)
                rowNames05{nameNum,1}=[rowNames05{nameNum,1},'_rank05'];
            end        
            tempResultsByNodeRank05=array2table(table2array(tempResultsByNodeRank05)','VariableNames',UseLabels,'RowNames',rowNames05);

            RefZ_RC_WholeCM=[tempResultsByCMRank50;tempResultsByCMRank25;tempResultsByCMRank10;tempResultsByCMRank05];
            RefZ_RC_WholeCM.Properties.RowNames={'Rank50','Rank25','Rank10','Rank05'};
            RefZ_RC_Degree=[tempResultsByNodeRank50;tempResultsByNodeRank25;tempResultsByNodeRank10;tempResultsByNodeRank05];
            RefZ_RC_Z=array2table(vertRSM2SymRSM(nan_zscore(mat2uppertriuvectormat(table2array(RefZ_RC))),[],1),'VariableNames',RefZ_RC.Properties.VariableNames,'RowNames',RefZ_RC.Properties.RowNames);
            RefZ_RF_Z = array2table(nan_zscore(table2array(RefZ_RF)')','VariableNames',RefZ_RF.Properties.VariableNames);       

            save(SaveNames{1,1},'Ref_RF','Ref_RC','RefZ_RF','RefZ_RC','Ref_RC_WholeCM','Ref_RC_Degree','Ref_RC_Z','Ref_RF_Z','RefZ_RC_WholeCM','RefZ_RC_Degree','RefZ_RC_Z','RefZ_RF_Z');
    
        end
    end
toc    
end
for parcelNum=1:numParcels
    ParcelName=ParcelNames{parcelNum,1};
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,ssAppend],'/',ParcelName,'/'];
    GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,ssAppend],'/',ParcelName,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,ssAppend],'/',ParcelName,'/'];
    
    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);     
    end
    if ~exist(GroupFiguresDir,'file')
        mkdir(GroupFiguresDir);     
    end
    if ~exist(GroupBrainMapsDir,'file')
        mkdir(GroupBrainMapsDir);     
    end     
    parcelName=ParcelName;
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
    if length(UseLabels)>1
        [RCMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','Ref_RC',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName);
        if istable(RCMats)
            RCMats=table2array(RCMats);
        end
        [RFMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','Ref_RF',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName); 
        if istable(RFMats)
            RFMats=table2array(RFMats);
        end    
        RCMatsVert=mat2uppertriuvectormat(RCMats);
        ssRSM_Ref_RF=single(corrcoef(RFMats));
        parcelRSM_Ref_RF=single(corrcoef(RFMats'));
        ssRSM_Ref_RF_scaled=single(corrcoef(scaleVals(RFMats,2)));
        parcelRSM_Ref_RF_scaled=single(corrcoef(scaleVals(RFMats',2)));  

        if length(UseLabels) <= 200
            ssRSM_Ref_RC=single(corrcoef(RCMatsVert));
            parcelRSM_Ref_RC=single(corrcoef(RCMatsVert'));
            ssRSM_Ref_RC_scaled=single(corrcoef(scaleVals(RCMatsVert,2)));
            parcelRSM_Ref_RC_scaled=single(corrcoef(scaleVals(RCMatsVert',2)));  
        else
            ssRSM_Ref_RC=[];
            parcelRSM_Ref_RC=[];
            ssRSM_Ref_RC_scaled=[];
            parcelRSM_Ref_RC_scaled=[];
        end

        [Group_Ref_RF,Group_Ref_RC,Group_Ref_RCMat_Ts,Group_Ref_RCMat_Means,Ref_RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(RFMats,RCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
            'vecName','Ref_RF','matName','Ref_RC');
        if ComputeReliability == 1
            Ref_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
            Ref_RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
            Ref_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
            Ref_RC_losoReliability=ReliabilityResults.Mat_losoReliability;
            save([GroupAnalysisDir,'Group_Ref_ParcelReliability'],'Ref_RF_SplitHalfReliability','Ref_RC_SplitHalfReliability','Ref_RF_losoReliability','Ref_RC_losoReliability'); 
        end
        save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
        save([GroupAnalysisDir,'Group_Ref_RF'],'Group_Ref_RF');
        save([GroupAnalysisDir,'Group_Ref_RC'],'Group_Ref_RC');
        save([GroupAnalysisDir,'ThirdOrderRSMs_Ref'],'ssRSM_Ref_RF','parcelRSM_Ref_RF','ssRSM_Ref_RF_scaled','parcelRSM_Ref_RF_scaled','ssRSM_Ref_RC','parcelRSM_Ref_RC','ssRSM_Ref_RC_scaled','parcelRSM_Ref_RC_scaled');
        save([GroupAnalysisDir,'Group_Ref_RCMat'],'Group_Ref_RCMat_Ts','Group_Ref_RCMat_Means','Ref_RC_Degree_Table');  

        [RCMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','RefZ_RC',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName);
        if istable(RCMats)
            RCMats=table2array(RCMats);
        end
        [RFMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','RefZ_RF',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName); 
        if istable(RFMats)
            RFMats=table2array(RFMats);
        end    
        RCMatsVert=mat2uppertriuvectormat(RCMats);
        ssRSM_RefZ_RF=single(corrcoef(RFMats));
        parcelRSM_RefZ_RF=single(corrcoef(RFMats'));
        ssRSM_RefZ_RF_scaled=single(corrcoef(scaleVals(RFMats,2)));
        parcelRSM_RefZ_RF_scaled=single(corrcoef(scaleVals(RFMats',2)));  

        if length(UseLabels) <= 200
            ssRSM_RefZ_RC=single(corrcoef(RCMatsVert));
            parcelRSM_RefZ_RC=single(corrcoef(RCMatsVert'));
            ssRSM_RefZ_RC_scaled=single(corrcoef(scaleVals(RCMatsVert,2)));
            parcelRSM_RefZ_RC_scaled=single(corrcoef(scaleVals(RCMatsVert',2)));  
        else
            ssRSM_RefZ_RC=[];
            parcelRSM_RefZ_RC=[];
            ssRSM_RefZ_RC_scaled=[];
            parcelRSM_RefZ_RC_scaled=[];
        end

        [Group_RefZ_RF,Group_RefZ_RC,Group_RefZ_RCMat_Ts,Group_RefZ_RCMat_Means,RefZ_RC_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(RFMats,RCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,'MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps,...
            'vecName','RefZ_RF','matName','RefZ_RC');
        if ComputeReliability == 1
            RefZ_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
            RefZ_RC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
            RefZ_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
            RefZ_RC_losoReliability=ReliabilityResults.Mat_losoReliability;
            save([GroupAnalysisDir,'Group_RefZ_ParcelReliability'],'RefZ_RF_SplitHalfReliability','RefZ_RC_SplitHalfReliability','RefZ_RF_losoReliability','RefZ_RC_losoReliability'); 
        end
        save([GroupAnalysisDir,'Group_RefZ_RF'],'Group_RefZ_RF');
        save([GroupAnalysisDir,'Group_RefZ_RC'],'Group_RefZ_RC');
        save([GroupAnalysisDir,'ThirdOrderRSMs_RefZ'],'ssRSM_RefZ_RF','parcelRSM_RefZ_RF','ssRSM_RefZ_RF_scaled','parcelRSM_RefZ_RF_scaled','ssRSM_RefZ_RC','parcelRSM_RefZ_RC','ssRSM_RefZ_RC_scaled','parcelRSM_RefZ_RC_scaled');
        save([GroupAnalysisDir,'Group_RefZ_RCMat'],'Group_RefZ_RCMat_Ts','Group_RefZ_RCMat_Means','RefZ_RC_Degree_Table');      
    else
        RCMats=[];
        [RFMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','Ref_RF',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName); 
        if istable(RFMats)
            RFMats=table2array(RFMats);
        end    
        RCMatsVert=[];
        ssRSM_Ref_RF=single(corrcoef(RFMats));
        parcelRSM_Ref_RF=single(corrcoef(RFMats'));
        ssRSM_Ref_RF_scaled=single(corrcoef(scaleVals(RFMats,2)));
        parcelRSM_Ref_RF_scaled=single(corrcoef(scaleVals(RFMats',2)));  

        ssRSM_Ref_RC=[];
        parcelRSM_Ref_RC=[];
        ssRSM_Ref_RC_scaled=[];
        parcelRSM_Ref_RC_scaled=[];
    
        [Group_Ref_RF,ReliabilityResults] = GroupVecParcellationSummaryFigures(RFMats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'vecName','Ref_RF');
        if ComputeReliability == 1
            Ref_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
            Ref_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
            save([GroupAnalysisDir,'Group_Ref_ParcelReliability'],'Ref_RF_SplitHalfReliability','Ref_RF_losoReliability'); 
        end   

        save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
        save([GroupAnalysisDir,'Group_Ref_RF'],'Group_Ref_RF');
        save([GroupAnalysisDir,'ThirdOrderRSMs_Ref'],'ssRSM_Ref_RF','parcelRSM_Ref_RF','ssRSM_Ref_RF_scaled','parcelRSM_Ref_RF_scaled');

        RCMats=[];
        [RFMats] = CompileND(ExperimentsDir,fmriprep_table,...
            'AnalysisType',AnalysisType,...
            'LoadVarName','RefZ_RF',...
            'DataType','ByParcellation',...
            'LoadVarFormat','Table',...
            'AnalysisName',AnalysisName,...
            'TableOrArray','Table',...
            'SubjectOrRun',SubjectOrRun,...
            'ParcelName',ParcelName); 
        if istable(RFMats)
            RFMats=table2array(RFMats);
        end    
        RCMatsVert=[];
        ssRSM_RefZ_RF=single(corrcoef(RFMats));
        parcelRSM_RefZ_RF=single(corrcoef(RFMats'));
        ssRSM_RefZ_RF_scaled=single(corrcoef(scaleVals(RFMats,2)));
        parcelRSM_RefZ_RF_scaled=single(corrcoef(scaleVals(RFMats',2)));  
        
        ssRSM_RefZ_RC=[];
        parcelRSM_RefZ_RC=[];
        ssRSM_RefZ_RC_scaled=[];
        parcelRSM_RefZ_RC_scaled=[];

        [Group_RefZ_RF,ReliabilityResults] = GroupVecParcellationSummaryFigures(RFMats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,'vecName','RefZ_RF');
        if ComputeReliability == 1
            RefZ_RF_SplitHalfReliability=ReliabilityResults.Vec_SplitHalfReliability;
            RefZ_RF_losoReliability=ReliabilityResults.Vec_losoReliability;
            save([GroupAnalysisDir,'Group_RefZ_ParcelReliability'],'RefZ_RF_SplitHalfReliability','RefZ_RF_losoReliability'); 
        end   

        save([GroupAnalysisDir,'AnalysisParameters'],'AnalysisParameters');
        save([GroupAnalysisDir,'Group_RefZ_RF'],'Group_RefZ_RF');
        save([GroupAnalysisDir,'ThirdOrderRSMs_RefZ'],'ssRSM_RefZ_RF','parcelRSM_RefZ_RF','ssRSM_RefZ_RF_scaled','parcelRSM_RefZ_RF_scaled');       
    end    
end

end

