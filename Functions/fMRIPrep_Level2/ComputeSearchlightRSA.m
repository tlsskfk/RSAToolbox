function [AnalysisParameters] = ComputeSearchlightRSA(fmriprep_table,ExperimentsDir,varargin)
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
[regressN] = VariableSetter('regressN',[],varargin);
[pRSM_FileName] = VariableSetter('pRSM_FileName',[],varargin);
[pRSM_VarNames] = VariableSetter('pRSM_VarNames',[],varargin);
%Subject or run level analysis. Will prompt request.
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);

%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);
if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            [filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType','RSA','TitleTextName','Select Analysis Parameters:');
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
    pRSM_FileName=AnalysisParameters.pRSM_FileName;
    pRSM_VarNames=AnalysisParameters.pRSM_VarNames; 
    FigureInfo=AnalysisParameters.FigureInfo;
    regressN=AnalysisParameters.regressN;
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


%% Compile filepaths for input files for the analysis
[filePaths_RSMs,~,RSMsName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RSMs','AnalysisName',RSMsName,'TitleTextName',['Select RSMs for',newline,'Representational Connectivity Analysis:']);
if isempty(ParcelNames)
    ParcelNames=filePaths_RSMs.Properties.VariableNames(:);
    ParcelNames=uiNameSelect(ParcelNames,'Select Searchlight to run.',1);
else
    ParcelNames=ParcelNames(ismember(ParcelNames,filePaths_RSMs.Properties.VariableNames(:)),:);
    ParcelNames=ParcelNames(ismember(ParcelNames,RefParcelNames),:);
end
if ~iscell(ParcelNames)
    ParcelNames={ParcelNames};
end
filePaths_RSMs=filePaths_RSMs(:,ParcelNames);
numParcels=length(ParcelNames);

if isempty(pRSM_FileName)
    [~,pRSM_FileNames]=getFolderAndFileNames('pRSMs/');
    SingleSelect=1;
    pRSM_FileName=uiNameSelect(pRSM_FileNames,'Select fmriprep_table used:',SingleSelect);
end
if iscell(pRSM_FileName)
    pRSM_FileName=pRSM_FileName{1,1};
end
load(['pRSMs/',pRSM_FileName],'pRSMs');
if isempty(pRSM_VarNames)
    SingleSelect=0;
    pRSM_VarNames=uiNameSelect(pRSMs.Properties.VariableNames(:),'Select fmriprep_table used:',SingleSelect);
end
pRSMs=table2array(pRSMs(:,pRSM_VarNames));

if isempty(regressN)
    regressN=uiEnterName('0','Enter # of regression analyses to run');
    regressN=str2num(regressN);
end

if regressN>0
    regGroups=cell(regressN,2);
    for regNum=1:regressN
        regGroups{regNum,1}=uiNameSelect(pRSM_VarNames,['Select pRSMs for Regression - ', num2str(regNum),' of ',num2str(regressN)]);
        regGroups{regNum,2}=pRSMs(:,ismember(pRSM_VarNames,regGroups{regNum,1}));
    end
else
    regGroups=[];
end
    
%Set analysis type and analysis name. These values will be used when saving
AnalysisType='RSA';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName([strrep(RSMsName,'RSMs_',[AnalysisType,'_'])],['Enter name for ',AnalysisType,newline,'analysis below:']);
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
AnalysisParameters.pRSM_FileName=pRSM_FileName;
AnalysisParameters.pRSM_VarNames=pRSM_VarNames;
AnalysisParameters.FigureInfo=FigureInfo;
AnalysisParameters.regressN=regressN;
%% BIDsTable loop: Iterate through BIDsTable and perform analysis
iniPercentComplete=0; %Used to display progress
save('Temp_ComputeParcelRSA_Params', 'AnalysisParameters');
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
    descript1='desc-RSA'; %set file description name
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
            TempLoadData = load(LoadPath_RSMs,'RSMs','rsm_mask');
            RSMs=single(TempLoadData.RSMs); 
            rsm_mask=TempLoadData.rsm_mask;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_RSMs]);
            continue
        end               
        %% Run analysis here!!       
        %load(['Parcellations/',ParcelName],'UseLabels')
        RSAMat=corr(pRSMs,RSMs,'type','spearman');  
        RSAMaps=repmat(rsm_mask,[1,1,1,size(pRSMs,2)])*0;
        for j= 1:size(pRSMs,2)
            tempMap=rsm_mask*0;
            tempMap(rsm_mask==1)=RSAMat(j,:);
            RSAMaps(:,:,:,j)=tempMap;
        end
        if ~isempty(regGroups)
            constant=ones(size(RSMs,1),1);
            regTs=cell(1,size(regGroups,1));
            regBs=cell(1,size(regGroups,1));         
            for regNum=1:size(regGroups,1)
                tempRegMapT=repmat(rsm_mask,[1,1,1,size(regGroups{regNum,2},2)])*0;
                tempRegMapB=repmat(rsm_mask,[1,1,1,size(regGroups{regNum,2},2)])*0;
                [tempTs,tempBs] = FastOLSRegress(RSMs,[regGroups{regNum,2},constant]);
                tempTs=tempTs(:,1:end-1);
                tempBs=tempBs(:,1:end-1);
                for j= 1:size(regGroups{regNum,2},2)
                    tempMapT=rsm_mask*0;
                    tempMapB=rsm_mask*0;
                    tempMapT(rsm_mask==1)=tempTs(:,j);
                    tempMapB(rsm_mask==1)=tempBs(:,j);
                    tempRegMapT(:,:,:,j)=tempMapT;
                    tempRegMapB(:,:,:,j)=tempMapB;
                end
                regTs{1,regNum}=tempRegMapT;
                regBs{1,regNum}=tempRegMapB;
            end
        else
            regGroups=[];
            regTs=[];
            regBs=[];
        end
        save(SaveNames{1,1},'RSAMaps','regGroups','regTs','regBs','pRSM_VarNames','AnalysisParameters');
    end     
toc    
end


for parcelNum=1:numParcels
    ParcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',ParcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
        GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',ParcelName,'/'];
    end
    if ~exist(GroupAnalysisDir,'file')
        mkdir(GroupAnalysisDir);     
    end
    if ~exist(GroupBrainMapsDir,'file')
        mkdir(GroupBrainMapsDir);     
    end     

    
    [RSAMaps] = CompileND(ExperimentsDir,fmriprep_table,...
        'AnalysisType',AnalysisType,...
        'LoadVarName','RSAMaps',...
        'DataType','ByParcellation',...
        'LoadVarFormat','Array',...
        'AnalysisName',AnalysisName,...
        'TableOrArray','Array',...
        'SubjectOrRun',SubjectOrRun,...
        'RowCompile','All',...
        'ParcelName',ParcelName);
    numSS=size(RSAMaps,5);
    RSAMaps(RSAMaps==0)=nan;
    ThreshMask=single(sum(single(~isnan(squeeze(RSAMaps(:,:,:,1,:)))),4)>round(0.75*numSS));
    ThreshMask=repmat(ThreshMask,[1,1,1,size(RSAMaps,4),size(RSAMaps,5)]);
    RSAMaps(ThreshMask==0)=nan;
    [groupRSA_t,groupRSA_z,groupRSA_mean]=getTval(atanh(RSAMaps),5);
    groupRSA_t(isinf(groupRSA_t))=nan;
    groupRSA_z(isinf(groupRSA_z))=nan;
    groupRSA_mean(isinf(groupRSA_mean))=nan;

    SaveBrik_3mmMNI(groupRSA_t,pRSM_VarNames,[GroupBrainMapsDir,ParcelName,'_RSA_SLMaps_T']);
    SaveBrik_3mmMNI(groupRSA_z,pRSM_VarNames,[GroupBrainMapsDir,ParcelName,'_RSA_SLMaps_Z']);
    SaveBrik_3mmMNI(groupRSA_mean,pRSM_VarNames,[GroupBrainMapsDir,ParcelName,'_RSA_SLMaps_mean']);
    save([GroupAnalysisDir,ParcelName,'_RSA_SLMaps'],'groupRSA_t','groupRSA_z','groupRSA_mean','AnalysisParameters');    
end

end


