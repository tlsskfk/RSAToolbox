function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = FC_IndvDiff(fmriprep_table,ExperimentsDir,varargin)
%Template function for data processing from the BIDsTable
%Overwrite previously saved files (default is no or 0; yes = 1)
if nargin==0
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = load_fmriprep_table;
else
    [fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
end

[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
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
AnalysisType='FC_IndvDiff';
%Allows you to set name for this particular analysis
FC_UseVar_Name=uiNameSelect({'FC','FC_fisherZ','FC_ZNorm','FC_fisherZ_ZNorm'},'Select FC measure to use.',1);
%% Compile filepaths for input files for the analysis

[filePaths_FC,~,FCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ParcelFC','TitleTextName',['Select FC values for',newline,'individual difference measures:']);
useIndicies(:,useIndicies>size(filePaths_FC,1))=[];
if isempty(AnalysisName)
    FC_IndvDiffName=strrep(FCName,'ParcelFC_','');
    AnalysisName=uiEnterName([FC_IndvDiffName,'_',FC_UseVar_Name],['Enter name for ',AnalysisType,newline,'analysis below:']);
end
[filePaths_beh_vars,~,~] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','beh_vars');
useIndicies(:,useIndicies>size(filePaths_beh_vars,1))=[];
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
ParcelNames=filePaths_FC.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'Searchlight'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_FC=filePaths_FC(:,ParcelNames);
numParcels=length(ParcelNames);

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.FCName=FCName;
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
    [~,labelPairs2,~]=labels2uppertriuvectorlabels(UseLabels);
    numFC=length(labelPairs2);
    %[ RC_ind,RC_reverseInd ] = LabelPair2Ind( UseLabels,labelPairs2,'_2_' );
    All_FC=[];
    for loadNum = useIndicies
        try
            load(filePaths_FC.(parcelName){loadNum,1},'FC');
            FC=array2table(double(table2array(FC(:,FC_UseVar_Name)))','VariableNames',labelPairs2);
        catch
            FC=array2table(nan(1,numFC,'double'),'VariableNames',labelPairs2);
            All_FC=[All_FC;FC];
            continue
        end
        try
            All_FC=[All_FC;FC];
        catch
            FC=array2table(nan(1,numFC,'double'),'VariableNames',labelPairs2);
            All_FC=[All_FC;FC];
        end
    end
    AnalysisParams.UseTable=fmriprep_table(useIndicies,:);
    AnalysisParams.parcelName=parcelName;
    AnalysisParams.lmem_formula=lmem_formula;
    AnalysisParams.Parcellation=load(['Parcellations/',parcelName]);
    AnalysisParams.FCPaths=strrepCell(filePaths_FC(useIndicies,:).(parcelName),ExperimentsDir,'');
    AnalysisParams.FCPaths=strrepCell(AnalysisParams.FCPaths,'\','/');
    %% Run mixed effects model 4 RF
    for behavVarNum=1:numvars
        disp([parcelName,': ',beh_var_names{behavVarNum,1}])
        beh_lmem_formula=strrep(lmem_formula,'behavVar',beh_var_names{behavVarNum,1});
        if ~any(~isnan(table2array(BaseTable(:,beh_var_names(behavVarNum,1)))))
            continue
        end
        
        % Run LMEM on RC data
        FC_Table=[BaseTable(:,[beh_var_names(behavVarNum,1);motion_var_names;RFX_Intercepts;RFX_Slopes])];
        Temp_Results_FC=cell(length(labelPairs2),1);
        Results_FC=[];
        parfor FCNum=1:length(labelPairs2)
            fc_lmem_formula=strrep(beh_lmem_formula,'brainVar',labelPairs2{FCNum,1});
            if ~any(~isnan(table2array(All_FC(:,labelPairs2{FCNum,1}))))
                continue
            end    
            try
                [Temp_Results_FC{FCNum,1}] = LME_PlusContinuousStats([FC_Table,All_FC(:,labelPairs2{FCNum,1})],fc_lmem_formula);  
            catch
                Temp_Results_FC{FCNum,1}=[];
            end    
        end 
        
        %Reformat data
        Results_FC_Confound.Intercept = [];
        for i = 1:length(motion_var_names)
            Results_FC_Confound.(motion_var_names{i,1})=[];
        end
        
        for FCNum=1:length(labelPairs2)
            if ~isempty(Temp_Results_FC{FCNum,1})                
                Results_FC=[Results_FC;Temp_Results_FC{FCNum,1}(labelPairs2{FCNum,1},:)];
                Results_FC_Confound.Intercept=[Results_FC_Confound.Intercept;Temp_Results_FC{FCNum,1}('Intercept',:)];
                Results_FC_Confound.Intercept.Properties.RowNames{FCNum,1}=labelPairs2{FCNum,1};
                for i = 1:length(motion_var_names)
                    Results_FC_Confound.(motion_var_names{i,1})=[Results_FC_Confound.(motion_var_names{i,1});Temp_Results_FC{FCNum,1}(motion_var_names{i,1},:)];
                    Results_FC_Confound.(motion_var_names{i,1}).Properties.RowNames{FCNum,1}=labelPairs2{FCNum,1};
                end
            else
                nanfill=array2table(nan(1,size(Results_FC,2)),'VariableNames',Results_FC.Properties.VariableNames);
                Results_FC=[Results_FC;nanfill];
                Results_FC.Properties.RowNames{FCNum,1}=labelPairs2{FCNum,1};
                Results_FC_Confound.Intercept=[Results_FC_Confound.Intercept;nanfill];
                Results_FC_Confound.Intercept.Properties.RowNames{FCNum,1}=labelPairs2{FCNum,1};
                for i = 1:length(motion_var_names)
                    Results_FC_Confound.(motion_var_names{i,1})=[Results_FC_Confound.(motion_var_names{i,1});nanfill];
                    Results_FC_Confound.(motion_var_names{i,1}).Properties.RowNames{FCNum,1}=labelPairs2{FCNum,1};
                end  
            end
        end
        if ~exist(GroupAnalysisDir,'file')
            mkdir(GroupAnalysisDir);
        end      
        save([GroupAnalysisDir,beh_var_names{behavVarNum,1}],'Results_FC','Results_FC_Confound');
        toc       
    end 
    save([GroupAnalysisDir,'AnalysisParams_IndvDiff'],'AnalysisParams');
end
end 

