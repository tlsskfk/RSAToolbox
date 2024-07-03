function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = RCA_SubSample_GroupAnalyses(fmriprep_table,ExperimentsDir,varargin)
%Written by David Rothlein
%Template function for data processing from the BIDsTable
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
%Subject or run level analysis. Will prompt request.
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
%Subject or run level analysis. Will prompt request.
[SSVarNames] = VariableSetter('SSVarNames',[],varargin);
[CMThresholds] = VariableSetter('CMThresholds',[0.001,0.005,0.01,0.05],varargin);
if isempty(fmriprep_table_name)
    [~,fmriprep_table_names]=getFolderAndFileNames('fmriprep_table/');
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect(fmriprep_table_names,'Select fmriprep_table used:',SingleSelect);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
end

%%Set up group Analysis directory
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'/Experiments/','');
GroupDir=[GroupDir,'/GroupAnalysis/'];
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');
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
AnalysisType='RCA_SubSample_Group';

%% Select / SubSample to use
[filePaths_subsamples,~,SubSampleName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','subsample','TitleTextName',['Select subsamples to use:']);
filePaths_subsamples=filePaths_subsamples.(SubSampleName);
AnalysisNameFilter=strrep(SubSampleName,'subsample_','');
LoadVarSubsample='SubAssignsByVol';
for dataInd=useIndicies
    try
        TempLoadDataSubsample = load(filePaths_subsamples{dataInd,1},LoadVarSubsample);
        numSubs=size(TempLoadDataSubsample.(LoadVarSubsample),2);
    catch
        continue
    end
    break
end
   
%% Compute SubSample motion Vars
[filePaths_CounfoundTCs] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','confounds','AnalysisName','ConfoundTC');
filePaths_CounfoundTCs=filePaths_CounfoundTCs.ConfoundTC;
LoadVarConfound='ConfoundTCs';
UseCleanVars={'framewise_displacement'};
numCleanVars=length(UseCleanVars);
CleanVarRemoveInd=zeros(length(useIndicies),1);
meanClean=nan(numSubs,numCleanVars,TotalRuns);
maxClean=nan(numSubs,numCleanVars,TotalRuns);
for dataInd=useIndicies
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    ConfoundData=cell(1);
    subsampleData=cell(1);   
    for run=1:numRuns
        loadInd=dataInd+run-1; %Set index to load        
        %try to load subsample data and skip if fails
        try
            TempLoadDataSubsample = load(filePaths_subsamples{loadInd,1},LoadVarSubsample);
            TempLoadDataSubsample = TempLoadDataSubsample.(LoadVarSubsample);
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVarSubsample,newline,filePaths_subsamples{loadInd,1}]);
            continue
        end             
        try
            TempLoadDataConfound = load(filePaths_CounfoundTCs{loadInd,1},LoadVarConfound);
            TempLoadDataConfound = TempLoadDataConfound.(LoadVarConfound);
            TempLoadDataConfound = TempLoadDataConfound(:,UseCleanVars);
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVarConfound,newline,filePaths_CounfoundTCs{loadInd,1}]);
            continue
        end
        ConfoundData{count,1}=TempLoadDataConfound;
        subsampleData{count,1}=TempLoadDataSubsample;
        count=count+1;
    end
    if count > 1
        [meanClean(:,:,dataInd),maxClean(:,:,dataInd)] = Compute_SubSampleVariable(ConfoundData,subsampleData); %(Subsample by Var)
    else
        CleanVarRemoveInd(dataInd,1)=1;
    end
end
clean_fmriprep_table=fmriprep_table;
clean_fmriprep_table(CleanVarRemoveInd==1,:)=[];
meanClean=permute(meanClean,[1,3,2]);
maxClean=permute(maxClean,[1,3,2]);
%% Select / identify parcellations to run
[filePaths_RCA,~,RCAAnalysisName] = BIDsDirSearch(ExperimentsDir,clean_fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','subsample_RCA','AnalysisNameFilter',AnalysisNameFilter,'TitleTextName',['Select RF and RC values for',newline,'individual difference measures:']);
ParcelNames=filePaths_RCA.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_RCA=filePaths_RCA(:,ParcelNames);
numParcels=length(ParcelNames);

RCAs2Run=uiNameSelect({'RF','RC'},'Select RCA values (select both for both)',0);

%% Compile SubSample Values

[subsample_data_all,clean_fmriprep_table,SubSampleInd] = CompileND(ExperimentsDir,clean_fmriprep_table,...
    'LoadVarName','subsample_data',...
    'LoadVarFormat','Array',...
    'DataType','ByParcellation',...
    'AnalysisType','subsample_RCA',...       
    'AnalysisName',RCAAnalysisName,...
    'SubjectOrRun',SubjectOrRun,...
    'TableOrArray','Array',...
    'ParcelName',ParcelNames{1,1});

meanClean(:,SubSampleInd==0,:)=[];
maxClean(:,SubSampleInd==0,:)=[];



%% Record analysis params

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.RCAName=RCAAnalysisName;
AnalysisParams.SubjectOrRun=SubjectOrRun;
AnalysisParams.fmriprep_table_name=fmriprep_table_name;

%% SelectVarType and motion correction type
IndVars=uiNameSelect({'SubSampleData','LinearStep'},'Select Y variable',0);
motion_var_names=uiNameSelect({[UseCleanVars{1,1},'_Mean'],[UseCleanVars{1,1},'_Max']},'Select motion variables to include',0);

%% Set up LME design and formula

RFX_Intercepts=uiNameSelect(clean_fmriprep_table.Properties.VariableNames,'Select additional random effect intercepts to include:',0);
RFX_Slopes=uiNameSelect(clean_fmriprep_table.Properties.VariableNames,'Select additional random effect slopes to include:',0);
AnalysisName=uiEnterName(RCAAnalysisName,['Enter name for ',AnalysisType,newline,'analysis below:']);
for IndVarNum = 1:size(IndVars,1)
    IndVar=IndVars{IndVarNum,1};
    lmem_formula=[IndVar,' ~ brainVar + '];
    if ~isempty(motion_var_names)
        for i = 1:length(motion_var_names)
            lmem_formula =  [lmem_formula,motion_var_names{i,1},' + '];
        end
    end
    if ~isempty(RFX_Intercepts)
        for i = 1:length(RFX_Intercepts)
            lmem_formula =  [lmem_formula,'(1|',RFX_Intercepts{i,1},')',' + '];
        end
    end  
    lmem_formula =  [lmem_formula,'(1 + brainVar|UniqueID) + '];
    if ~isempty(RFX_Slopes)
        for i = 1:length(RFX_Slopes)
            lmem_formula =  [lmem_formula,'(1 + brainVar|',RFX_Slopes{i,1},')',' + '];
        end
    end  

    if strcmp(lmem_formula(1,end-2:end),' + ')
        lmem_formula(:,end-2:end)=[];
    end

    %BaseTable=[fmriprep_table(useIndicies,RFX_Intercepts),fmriprep_table(useIndicies,RFX_Slopes),All_MotionStats,All_Variables];
    RFMaps=cell(numParcels,1);
    RCMaps=cell(numParcels,5,length(CMThresholds));
    CircPlots=cell(numParcels,1);
    for parcelNum=1:numParcels
        Results_RC=[];
        Results_RF=[];
        CircData=[];
        parcelName=ParcelNames{parcelNum,1};
        if bySS == 1
            GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
            GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
            GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
        else
            GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
            GroupFiguresDir=[FigureDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
            GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
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
        
        load(['Parcellations/',parcelName],'UseMask','UseLabels');
        UseLabels=UseLabels(:);
        disp([IndVar,': ',parcelName])
        
        if any(ismember(RCAs2Run,'RF'))
            tic
            disp('Initializing RF')
            numRF=length(UseLabels);
            AnalysisParams.parcelName=parcelName;
            AnalysisParams.lmem_formula=lmem_formula;
            AnalysisParams.Parcellation=load(['Parcellations/',parcelName]);

            %% Compile Subsample 
            [RCAData,compiled_fmriprep_table,selectind_RF] = CompileND(ExperimentsDir,clean_fmriprep_table,...
                'LoadVarName','subRCA_RF',...
                'LoadVarFormat','Array',...
                'DataType','ByParcellation',...
                'AnalysisType','subsample_RCA',...       
                'AnalysisName',RCAAnalysisName,...
                'SubjectOrRun',SubjectOrRun,...
                'TableOrArray','Array',...
                'ParcelName',parcelName); 
            RCAData=permute(RCAData,[1,3,2]); 
            RCADataTable=[];
            for RCADatai = 1:size(RCAData,3)
                tempRCAData=RCAData(:,:,RCADatai);
                RCADataTable=[RCADataTable,array2table(double(tempRCAData(:)),'VariableNames',UseLabels(RCADatai))];
            end

            %% Set up basetable for RF
            UniqueID=repmat(strrepCell(join([compiled_fmriprep_table.experiment,compiled_fmriprep_table.sub]),' ','_')',[numSubs,1]); 
            UniqueID=cell2table(UniqueID(:),'VariableNames',{'UniqueID'});
            Use_meanClean=meanClean(:,selectind_RF==1,:);
            Use_meanClean=array2table(double(Use_meanClean(:)),'VariableNames',{[UseCleanVars{1,1},'_Mean']});
            Use_maxClean=maxClean(:,selectind_RF==1,:);  
            Use_maxClean=array2table(double(Use_maxClean(:)),'VariableNames',{[UseCleanVars{1,1},'_Max']});
            intTable=[];
            slpTable=[];
            if ~isempty(RFX_Intercepts)            
                for iInt = 1:length(RFX_Intercepts)
                    tempVar=repmat(compiled_fmriprep_table.(RFX_Intercepts{iInt,1})',[numSubs,1]);
                    intTable=[intTable,cell2table(tempVar(:),'VariableNames',RFX_Intercepts(iInt,1))];
                end
            end  

            if ~isempty(RFX_Slopes)
                for iSlp = 1:length(RFX_Slopes)
                    tempVar=repmat(compiled_fmriprep_table.(RFX_Slopes{iSlp,1})',[numSubs,1]);
                    slpTable=[slpTable,cell2table(tempVar(:),'VariableNames',RFX_Slopes(iSlp,1))];
                end
            end   
            if strcmpi(IndVar,'SubSampleData')
                use_subsample_data=subsample_data_all(:,selectind_RF==1);
                indVarTable=array2table(double(use_subsample_data(:)),'VariableNames',{IndVar});
            else
                stepData=repmat([1:numSubs]',[1,height(compiled_fmriprep_table)]);
                indVarTable=array2table(double(stepData(:)),'VariableNames',{IndVar});
            end
            BaseTable=[UniqueID,indVarTable,Use_meanClean,Use_maxClean,intTable,slpTable];
            Temp_Results_RF=cell(length(UseLabels),1);
            %% Run LMEM on RF data
            toc
            tic
            disp('Running RF LME')
            parfor_progress(length(UseLabels));
            parfor RFNum=1:length(UseLabels)
                rf_lmem_formula=strrep(lmem_formula,'brainVar',UseLabels{RFNum,1});
                if ~any(~isnan(table2array(RCADataTable(:,UseLabels{RFNum,1}))))
                    continue
                end 
                try
                    [Temp_Results_RF{RFNum,1}] = LME_PlusContinuousStats([BaseTable,RCADataTable(:,UseLabels{RFNum,1})],rf_lmem_formula);
                catch
                    Temp_Results_RF{RFNum,1}=[];
                end
                parfor_progress;
            end 
            parfor_progress(0);
            toc
            tic
            disp('Compiling RF Results')       
            Results_RF_Confound.Intercept = [];
            for i = 1:length(motion_var_names)
                Results_RF_Confound.(motion_var_names{i,1})=[];
            end   
            Results_RF=[];
            for RFNum=1:length(UseLabels)
                if ~isempty(Temp_Results_RF{RFNum,1})
                    Results_RF=[Results_RF;Temp_Results_RF{RFNum,1}(UseLabels{RFNum,1},:)];
                    Results_RF_Confound.Intercept=[Results_RF_Confound.Intercept;Temp_Results_RF{RFNum,1}('Intercept',:)];
                    Results_RF_Confound.Intercept.Properties.RowNames{RFNum,1}=UseLabels{RFNum,1};
                    for i = 1:length(motion_var_names)
                        Results_RF_Confound.(motion_var_names{i,1})=[Results_RF_Confound.(motion_var_names{i,1});Temp_Results_RF{RFNum,1}(motion_var_names{i,1},:)];
                        Results_RF_Confound.(motion_var_names{i,1}).Properties.RowNames{RFNum,1}=UseLabels{RFNum,1};
                    end
                else
                    nanfill=array2table(nan(1,size(Results_RF,2)),'VariableNames',Results_RF.Properties.VariableNames);
                    Results_RF=[Results_RF;nanfill];
                    Results_RF.Properties.RowNames{RFNum,1}=UseLabels{RFNum,1};
                    Results_RF_Confound.Intercept=[Results_RF_Confound.Intercept;nanfill];
                    Results_RF_Confound.Intercept.Properties.RowNames{RFNum,1}=UseLabels{RFNum,1};
                    for i = 1:length(motion_var_names)
                        Results_RF_Confound.(motion_var_names{i,1})=[Results_RF_Confound.(motion_var_names{i,1});nanfill];
                        Results_RF_Confound.(motion_var_names{i,1}).Properties.RowNames{RFNum,1}=UseLabels{RFNum,1};
                    end  
                end
            end  
            mapName=['RF_',IndVar,'_',parcelName];
            mapVals=p2z(Results_RF.lme_pValue,2).*sign(Results_RF.lme_tStat);
            [ RFMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(RFMaps,{mapName},[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'RFMaps');
            
            if size(mapVals,1)<50
                CircData.rVec=Results_RF.lme_Estimate;
                CircData.pVec=Results_RF.lme_pValue;
                if size(mapVals,1) < 15
                    CircData.FontSize = 16;
                elseif size(mapVals,1) >= 15 && size(mapVals,1) < 25
                    CircData.FontSize = 14;
                else
                    CircData.FontSize = 12;
                end
            else
                CircData=[];
            end
            toc
        end   
        if any(ismember(RCAs2Run,'RC'))
        %% Run 
            tic
            disp('Initializing RC')    
            [~,labelPairs]=labels2uppertriuvectorlabels(UseLabels);
            numRC=length(labelPairs);    
            AnalysisParams.parcelName=parcelName;
            AnalysisParams.lmem_formula=lmem_formula;
            AnalysisParams.Parcellation=load(['Parcellations/',parcelName]);
            %% Compile Subsample 
            [RCAData,compiled_fmriprep_table,selectind_RC] = CompileND(ExperimentsDir,clean_fmriprep_table,...
                'LoadVarName','subRCA_RC',...
                'LoadVarFormat','Array',...
                'DataType','ByParcellation',...
                'AnalysisType','subsample_RCA',...       
                'AnalysisName',RCAAnalysisName,...
                'SubjectOrRun',SubjectOrRun,...
                'TableOrArray','Array',...
                'ParcelName',parcelName); 
            RCAData=permute(RCAData,[1,3,2]); 
            RCADataTable=[];
            for RCADatai = 1:size(RCAData,3)
                tempRCAData=RCAData(:,:,RCADatai);
                RCADataTable=[RCADataTable,array2table(double(tempRCAData(:)),'VariableNames',labelPairs(RCADatai))];
            end

            %% Set up basetable for RF
            UniqueID=repmat(strrepCell(join([compiled_fmriprep_table.experiment,compiled_fmriprep_table.sub]),' ','_')',[numSubs,1]); 
            UniqueID=cell2table(UniqueID(:),'VariableNames',{'UniqueID'});
            Use_meanClean=meanClean(:,selectind_RC==1,:);
            Use_meanClean=array2table(double(Use_meanClean(:)),'VariableNames',{[UseCleanVars{1,1},'_Mean']});
            Use_maxClean=maxClean(:,selectind_RC==1,:);  
            Use_maxClean=array2table(double(Use_maxClean(:)),'VariableNames',{[UseCleanVars{1,1},'_Max']});
            intTable=[];
            slpTable=[];
            if ~isempty(RFX_Intercepts)            
                for iInt = 1:length(RFX_Intercepts)
                    tempVar=repmat(compiled_fmriprep_table.(RFX_Intercepts{iInt,1})',[numSubs,1]);
                    intTable=[intTable,cell2table(tempVar(:),'VariableNames',RFX_Intercepts(iInt,1))];
                end
            end  

            if ~isempty(RFX_Slopes)
                for iSlp = 1:length(RFX_Slopes)
                    tempVar=repmat(compiled_fmriprep_table.(RFX_Slopes{iSlp,1})',[numSubs,1]);
                    slpTable=[slpTable,cell2table(tempVar(:),'VariableNames',RFX_Slopes(iSlp,1))];
                end
            end   
            if strcmpi(IndVar,'SubSampleData')
                use_subsample_data=subsample_data_all(:,selectind_RC==1);
                indVarTable=array2table(double(use_subsample_data(:)),'VariableNames',{IndVar});
            else
                stepData=repmat([1:numSubs]',[1,height(compiled_fmriprep_table)]);
                indVarTable=array2table(double(stepData(:)),'VariableNames',{IndVar});
            end
            BaseTable=[UniqueID,indVarTable,Use_meanClean,Use_maxClean,intTable,slpTable];
            Temp_Results_RC=cell(length(labelPairs),1);
            toc
            %% Run LMEM on RC data
            tic
            disp('Running RC LME')
            parfor_progress(length(labelPairs));        
            parfor RCNum=1:length(labelPairs)
                rc_lmem_formula=strrep(lmem_formula,'brainVar',labelPairs{RCNum,1});
                if ~any(~isnan(table2array(RCADataTable(:,labelPairs{RCNum,1}))))
                    continue
                end 
                try
                    [Temp_Results_RC{RCNum,1}] = LME_PlusContinuousStats([BaseTable,RCADataTable(:,labelPairs{RCNum,1})],rc_lmem_formula);
                catch
                    Temp_Results_RC{RCNum,1}=[];
                end
                parfor_progress;
            end 
            parfor_progress(0);
            toc
            tic
            disp('Compiling RC Results')         
            Results_RC_Confound.Intercept = [];
            for i = 1:length(motion_var_names)
                Results_RC_Confound.(motion_var_names{i,1})=[];
            end   

            for RCNum=1:length(labelPairs)
                if ~isempty(Temp_Results_RC{RCNum,1})                
                    Results_RC=[Results_RC;Temp_Results_RC{RCNum,1}(labelPairs{RCNum,1},:)];
                    Results_RC_Confound.Intercept=[Results_RC_Confound.Intercept;Temp_Results_RC{RCNum,1}('Intercept',:)];
                    Results_RC_Confound.Intercept.Properties.RowNames{RCNum,1}=labelPairs{RCNum,1};
                    for i = 1:length(motion_var_names)
                        Results_RC_Confound.(motion_var_names{i,1})=[Results_RC_Confound.(motion_var_names{i,1});Temp_Results_RC{RCNum,1}(motion_var_names{i,1},:)];
                        Results_RC_Confound.(motion_var_names{i,1}).Properties.RowNames{RCNum,1}=labelPairs{RCNum,1};
                    end
                else
                    nanfill=array2table(nan(1,size(Results_RC,2)),'VariableNames',Results_RC.Properties.VariableNames);
                    Results_RC=[Results_RC;nanfill];
                    Results_RC.Properties.RowNames{RCNum,1}=labelPairs{RCNum,1};
                    Results_RC_Confound.Intercept=[Results_RC_Confound.Intercept;nanfill];
                    Results_RC_Confound.Intercept.Properties.RowNames{RCNum,1}=labelPairs{RCNum,1};
                    for i = 1:length(motion_var_names)
                        Results_RC_Confound.(motion_var_names{i,1})=[Results_RC_Confound.(motion_var_names{i,1});nanfill];
                        Results_RC_Confound.(motion_var_names{i,1}).Properties.RowNames{RCNum,1}=labelPairs{RCNum,1};
                    end  
                end
            end
            [ pMat ] = vertRSM2SymRSM( Results_RC.lme_pValue );
            [ rMat ] = vertRSM2SymRSM( Results_RC.lme_Estimate );
            MapLabels=[];
            MapVals=[];
            for pThresh=1:length(CMThresholds)
                [Results] = ConnectivityMatrixStats(rMat,'pmat',pMat,'Thresholds',CMThresholds(1,pThresh));
                tempVarNames=Results.Properties.VariableNames(:)';
                for vNum = 1:length(tempVarNames)
                    tempVarNames{1,vNum}=[tempVarNames{1,vNum},'_',num2str4filename(CMThresholds(1,pThresh),4)];
                end
                MapLabels=[MapLabels,tempVarNames];
                MapVals=[MapVals,table2array(Results)];
            end
            [ RCMaps ] = MakeParcellationMap( MapVals,UseMask);
            mapName=['RC_',IndVar,'_',parcelName];
            SaveBrik_3mmMNI(RCMaps,MapLabels,[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'RCMaps');            
            RC_Stats_Table=array2table(MapVals,'VariableNames',MapLabels);
            RC_Stats_Table.Properties.RowNames=UseLabels;
            
            if size(mapVals,1)<50
                CircData.rMat=rMat;
                CircData.pMat=pMat;
                if size(mapVals,1) < 15
                    CircData.FontSize = 12;
                elseif size(mapVals,1) >= 15 && size(mapVals,1) < 25
                    CircData.FontSize = 10;
                else
                    CircData.FontSize = 8;
                end
            else
                CircData=[];
            end            
            
            if ~isempty(CircData) 
                if isfield(CircData,'rMat') && isfield(CircData,'rVec')
                    figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(CircData.rVec,CircData.pVec,CircData.rMat,CircData.pMat,strrepCell(UseLabels,'_',' '),'FontSize',CircData.FontSize,'pThresh',0.05,'pThreshVec',0.05);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,IndVar,'_',parcelName,'p05.png'],'-Transparent','-png','-m5')
                    end
                    close   
                    
                    figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(CircData.rVec,CircData.pVec,CircData.rMat,CircData.pMat,strrepCell(UseLabels,'_',' '),'FontSize',CircData.FontSize,'pThresh',0.005,'pThreshVec',0.005);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,IndVar,'_',parcelName,'p005.png'],'-Transparent','-png','-m5')
                    end    
                    close
                    
                    figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(CircData.rVec,CircData.pVec,CircData.rMat,CircData.pMat,strrepCell(UseLabels,'_',' '),'FontSize',CircData.FontSize);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,IndVar,'_',parcelName,'_bonf.png'],'-Transparent','-png','-m5')
                    end
                    close                    
                end
            end
            toc
        end 
        disp('Saving Results')
        tic

        SaveName=[GroupAnalysisDir,IndVar,'.mat'];
        save([GroupAnalysisDir,'AnalysisParams_SubSample'],'AnalysisParams');
        if ~exist(SaveName,'file')
            if ~isempty(Results_RF) && ~isempty(Results_RC)
                save(SaveName,'Results_RF','Results_RF_Confound','Results_RC','Results_RC_Confound','RC_Stats_Table');
            elseif ~isempty(Results_RF)
                save(SaveName,'Results_RF','Results_RF_Confound','RC_Stats_Table');
            elseif ~isempty(Results_RC)
                save(SaveName,'Results_RC','Results_RC_Confound','RC_Stats_Table');  
            end
        else
            if ~isempty(Results_RF) && ~isempty(Results_RC)
                save(SaveName,'Results_RF','Results_RF_Confound','Results_RC','Results_RC_Confound');
            elseif ~isempty(Results_RF)
                if any(ismember(who('-file', SaveName),'Results_RC'))
                    load(SaveName,'Results_RC','Results_RC_Confound','RC_Stats_Table')
                    save(SaveName,'Results_RF','Results_RF_Confound','Results_RC','Results_RC_Confound','RC_Stats_Table');
                else
                    save(SaveName,'Results_RF','Results_RF_Confound');
                end
            elseif ~isempty(Results_RC)
                if any(ismember(who('-file', SaveName),'Results_RF'))
                    load(SaveName,'Results_RF','Results_RF_Confound')
                    save(SaveName,'Results_RF','Results_RF_Confound','Results_RC','Results_RC_Confound','RC_Stats_Table');
                else
                    save(SaveName,'Results_RC','Results_RC_Confound','RC_Stats_Table');
                end
            end    
        end
        toc       
    end  
end
end
 

