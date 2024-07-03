function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = RCA_IndvDiff(fmriprep_table,ExperimentsDir,varargin)
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

if isempty(fmriprep_table_name)
    [~,fmriprep_table_names]=getFolderAndFileNames('fmriprep_table/');
    SingleSelect=1;
    fmriprep_table_name=uiNameSelect(fmriprep_table_names,'Select fmriprep_table used:',SingleSelect);
    fmriprep_table_name=strrep(fmriprep_table_name,'.mat','');
end
ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'/Experiments/','');
GroupDir=[GroupDir,'/GroupAnalysis/'];
GroupDir=strrep(GroupDir,'//','/');
AnalysisParams=struct;

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
AnalysisType='RCA_IndvDiff';
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


%% Compile filepaths for input files for the analysis
[filePaths_RCA,~,RCAName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','RCA','TitleTextName',['Select RF and RC values for',newline,'individual difference measures:']);
if isempty(AnalysisName)
    RCA_IndvDiffName=strrep(RCAName,['RCA_',fmriprep_table_name,'_'],'');
    AnalysisName=uiEnterName(RCA_IndvDiffName,['Enter name for ',AnalysisType,newline,'analysis below:']);
end
[filePaths_beh_vars,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','beh_vars');
beh_var_names=[];
for i = useIndicies
    try
        load(filePaths_beh_vars.Overall{i,1});
    catch
        continue
    end
    beh_var_names=unique([beh_var_names;beh_vars.Properties.VariableNames(:)]);
end
beh_var_names=uiNameSelect(beh_var_names,'Select behavioral varaibles to compare.',0);
if ~iscell(beh_var_names)
    beh_var_names={beh_var_names};
end
if isempty(SSVarNames)
    SingleSelect=1; %Allows only a single value to be selected.
    [IncludeSSVars] = uiNameSelect({'Yes','No'},'Include SS vars?',SingleSelect);
    if strcmpi(IncludeSSVars,'Yes')
        IncludeSSVars=1;
    else
        IncludeSSVars=0;
    end
else
    [filePaths_ss_vars,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SSVars');
    IncludeSSVars=0;
end
if IncludeSSVars==1
    [filePaths_ss_vars,~,SSVarTitle] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','SSVars');
    for i = useIndicies
        try
            load(filePaths_ss_vars.(SSVarTitle){i,1},'SS_Vars');
        catch
            continue
        end
        SSVarNames=unique([SSVarNames;SS_Vars.Properties.VariableNames(:)]);
    end  
    SSVarNames=uiNameSelect(SSVarNames,'Select SS varaibles to compare.',0);
end

numvars=size(beh_var_names,1)+size(SSVarNames,1);

%% Identify parcellations availible and select ones to run analysis on.
ParcelNames=filePaths_RCA.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_RCA=filePaths_RCA(:,ParcelNames);
numParcels=length(ParcelNames);

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.RCAName=RCAName;
AnalysisParams.SubjectOrRun=SubjectOrRun;
AnalysisParams.fmriprep_table_name=fmriprep_table_name;


%% Load Behavioral variables
All_Variables=[];
if isempty(SSVarNames)
    for i = useIndicies
        try
            load(filePaths_beh_vars.Overall{i,1},'beh_vars');
        catch
            beh_vars=array2table(nan(1,numvars),'VariableNames',beh_var_names);
            All_Variables=[All_Variables;beh_vars];
            continue
        end
        try
            All_Variables=[All_Variables;beh_vars(:,beh_var_names)];
        catch
            beh_vars=array2table(nan(1,numvars),'VariableNames',beh_var_names);
            All_Variables=[All_Variables;beh_vars];
        end
    end
else
    if ~iscell(SSVarNames)
        SSVarNames={SSVarNames};
    end
    for i = useIndicies
        try
            load(filePaths_beh_vars.Overall{i,1},'beh_vars');
            beh_vars=beh_vars(:,beh_var_names);
        catch
            beh_vars=array2table(nan(1,length(beh_var_names)),'VariableNames',beh_var_names);
        end
        try
            load(filePaths_ss_vars.(SSVarTitle){i,1},'SS_Vars');
            SS_Vars=SS_Vars(:,SSVarNames);
        catch
            SS_Vars=array2table(nan(1,length(SSVarNames)),'VariableNames',SSVarNames);
        end       
        
        try
            All_Variables=[All_Variables;[beh_vars,SS_Vars]];
        catch
            beh_vars=array2table(nan(1,numvars),'VariableNames',beh_var_names);
            All_Variables=[All_Variables;array2table(nan(1,numvars),'VariableNames',[beh_var_names(:)',SSVarNames(:)'])];
        end
    end
end
All_Variables=array2table(double(table2array(All_Variables)),'VariableNames',All_Variables.Properties.VariableNames);
beh_var_names=All_Variables.Properties.VariableNames(:);
numvars=size(beh_var_names,1);


%% Select motion or confound variables for mixed effects model
if bySS==0
    motion_var_names=uiNameSelect({'MotionByRun_rot_mean','MotionByRun_rot_max','MotionByRun_trans_mean','MotionByRun_trans_max','MotionByRun_nsso_sum'},'Select motion variables to include',0);
else
    motion_var_names=uiNameSelect({'MotionBySS_rot_mean','MotionBySS_rot_max','MotionBySS_trans_mean','MotionBySS_trans_max','MotionBySS_nsso_sum'},'Select motion variables to include',0);
end
All_MotionStats=fmriprep_table(useIndicies,motion_var_names);
All_MotionStats=array2table(double(table2array(All_MotionStats)),'VariableNames',All_MotionStats.Properties.VariableNames);
RFX_Intercepts=uiNameSelect(fmriprep_table.Properties.VariableNames,'Select random effect intercepts to include:',0);
RFX_Slopes=uiNameSelect(fmriprep_table.Properties.VariableNames,'Select random effect slopes to include:',0);
lmem_formula=['behavVar ~ brainVar + '];
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
if ~isempty(RFX_Slopes)
    for i = 1:length(RFX_Slopes)
        lmem_formula =  [lmem_formula,'(1 + brainVar|',RFX_Slopes{i,1},')',' + '];
    end
end  

if strcmp(lmem_formula(1,end-2:end),' + ')
    lmem_formula(:,end-2:end)=[];
end
BaseTable=[fmriprep_table(useIndicies,RFX_Intercepts),fmriprep_table(useIndicies,RFX_Slopes),All_MotionStats,All_Variables];
for parcelNum=1:numParcels
    tic
    parcelName=ParcelNames{parcelNum,1};
    if bySS == 1
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',[AnalysisName,'_bySS'],'/',parcelName,'/'];
    else
        GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',parcelName,'/'];
    end
    load(['Parcellations/',parcelName],'UseLabels');
    numRF=length(UseLabels);
    [labelPairs1,~,~]=labels2uppertriuvectorlabels(UseLabels);
    numRC=length(labelPairs1);
    %[ RC_ind,RC_reverseInd ] = LabelPair2Ind( UseLabels,labelPairs1,'_2_' );
    
    All_RCA_RF=[];
    All_RCA_RC=[];
    for loadNum = useIndicies
        try
            load(filePaths_RCA.(parcelName){loadNum,1},'RCA_RF','RCA_RC');
            RCA_RF=array2table(double(table2array(RCA_RF)),'VariableNames',RCA_RF.Properties.VariableNames);
            RCA_RC=array2table(mat2uppertriuvectormat(double(table2array(RCA_RC)))','VariableNames',labelPairs1);
        catch
            RCA_RF=array2table(nan(1,numRF,'double'),'VariableNames',All_RCA_RF.Properties.VariableNames);
            RCA_RC=array2table(nan(1,numRC,'double'),'VariableNames',labelPairs1);
            All_RCA_RF=[All_RCA_RF;RCA_RF];
            All_RCA_RC=[All_RCA_RC;RCA_RC];
            continue
        end
        try
            All_RCA_RF=[All_RCA_RF;RCA_RF];
            All_RCA_RC=[All_RCA_RC;RCA_RC];
        catch
            RCA_RF=array2table(nan(1,numRF,'double'),'VariableNames',All_RCA_RF.Properties.VariableNames);
            RCA_RC=array2table(nan(1,numRC,'double'),'VariableNames',labelPairs1);
            All_RCA_RF=[All_RCA_RF;RCA_RF];
            All_RCA_RC=[All_RCA_RC;RCA_RC];
        end
    end
    AnalysisParams.UseTable=fmriprep_table(useIndicies,:);
    AnalysisParams.parcelName=parcelName;
    AnalysisParams.lmem_formula=lmem_formula;
    AnalysisParams.Parcellation=load(['Parcellations/',parcelName]);
    AnalysisParams.RCAPaths=strrepCell(filePaths_RCA(useIndicies,:).(parcelName),ExperimentsDir,'');
    AnalysisParams.RCAPaths=strrepCell(AnalysisParams.RCAPaths,'\','/');
    %% Run mixed effects model 4 RF
    for behavVarNum=1:numvars
        disp([parcelName,': ',beh_var_names{behavVarNum,1}])
        beh_lmem_formula=strrep(lmem_formula,'behavVar',beh_var_names{behavVarNum,1});
        RF_Table=[BaseTable(:,[beh_var_names(behavVarNum,1);motion_var_names;RFX_Intercepts;RFX_Slopes])];
        Temp_Results_RF=cell(length(UseLabels),1);
        Results_RF=[];
        if ~any(~isnan(table2array(BaseTable(:,beh_var_names(behavVarNum,1)))))
            continue
        end
            
        % Run LMEM on RF data
        parfor RFNum=1:length(UseLabels)
            rf_lmem_formula=strrep(beh_lmem_formula,'brainVar',UseLabels{RFNum,1});
            if ~any(~isnan(table2array(All_RCA_RF(:,UseLabels{RFNum,1}))))
                continue
            end 
            try
                [Temp_Results_RF{RFNum,1}] = LME_PlusContinuousStats([RF_Table,All_RCA_RF(:,UseLabels{RFNum,1})],rf_lmem_formula);
            catch
                Temp_Results_RF{RFNum,1}=[];
            end
        end
        
        % Run LMEM on RC data
        RC_Table=[BaseTable(:,[beh_var_names(behavVarNum,1);motion_var_names;RFX_Intercepts;RFX_Slopes])];
        Temp_Results_RC=cell(length(labelPairs1),1);
        Results_RC=[];
        parfor RCNum=1:length(labelPairs1)
            rc_lmem_formula=strrep(beh_lmem_formula,'brainVar',labelPairs1{RCNum,1});
            if ~any(~isnan(table2array(All_RCA_RC(:,labelPairs1{RCNum,1}))))
                continue
            end    
            try
                [Temp_Results_RC{RCNum,1}] = LME_PlusContinuousStats([RC_Table,All_RCA_RC(:,labelPairs1{RCNum,1})],rc_lmem_formula);  
            catch
                Temp_Results_RC{RCNum,1}=[];
            end    
        end 
        
        %Reformat data
        Results_RF_Confound.Intercept = [];
        for i = 1:length(motion_var_names)
            Results_RF_Confound.(motion_var_names{i,1})=[];
        end
        Results_RC_Confound=Results_RF_Confound;  
        
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
        for RCNum=1:length(labelPairs1)
            if ~isempty(Temp_Results_RC{RCNum,1})                
                Results_RC=[Results_RC;Temp_Results_RC{RCNum,1}(labelPairs1{RCNum,1},:)];
                Results_RC_Confound.Intercept=[Results_RC_Confound.Intercept;Temp_Results_RC{RCNum,1}('Intercept',:)];
                Results_RC_Confound.Intercept.Properties.RowNames{RCNum,1}=labelPairs1{RCNum,1};
                for i = 1:length(motion_var_names)
                    Results_RC_Confound.(motion_var_names{i,1})=[Results_RC_Confound.(motion_var_names{i,1});Temp_Results_RC{RCNum,1}(motion_var_names{i,1},:)];
                    Results_RC_Confound.(motion_var_names{i,1}).Properties.RowNames{RCNum,1}=labelPairs1{RCNum,1};
                end
            else
                nanfill=array2table(nan(1,size(Results_RC,2)),'VariableNames',Results_RC.Properties.VariableNames);
                Results_RC=[Results_RC;nanfill];
                Results_RC.Properties.RowNames{RCNum,1}=labelPairs1{RCNum,1};
                Results_RC_Confound.Intercept=[Results_RC_Confound.Intercept;nanfill];
                Results_RC_Confound.Intercept.Properties.RowNames{RCNum,1}=labelPairs1{RCNum,1};
                for i = 1:length(motion_var_names)
                    Results_RC_Confound.(motion_var_names{i,1})=[Results_RC_Confound.(motion_var_names{i,1});nanfill];
                    Results_RC_Confound.(motion_var_names{i,1}).Properties.RowNames{RCNum,1}=labelPairs1{RCNum,1};
                end  
            end
        end
        if ~exist(GroupAnalysisDir,'file')
            mkdir(GroupAnalysisDir);
        end
        save([GroupAnalysisDir,'AnalysisParams_IndvDiff'],'AnalysisParams');
        save([GroupAnalysisDir,beh_var_names{behavVarNum,1}],'Results_RF','Results_RF_Confound','Results_RC','Results_RC_Confound');
        toc       
    end  
end
end 

