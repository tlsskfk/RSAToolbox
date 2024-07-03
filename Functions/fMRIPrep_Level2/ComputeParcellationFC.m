function [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeParcellationFC(fmriprep_table,ExperimentsDir,varargin)
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
[FigureInfo] = VariableSetter('FigureInfo',[],varargin);
if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end

ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
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
[filePaths_ParcelTC,~,ParcelTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun',SubjectOrRun,'AnalysisType','ParcelTC','TitleTextName','Select timecourses for functional connenctome:');
ParcelNames=filePaths_ParcelTC.Properties.VariableNames;
ParcelNames=ParcelNames(:);
ParcelNames(ismember(ParcelNames,'WholeBrain'),:)=[];
ParcelNames=uiNameSelect(ParcelNames,'Select parcellations to run.');
filePaths_ParcelTC=filePaths_ParcelTC(:,ParcelNames);
numParcels=length(ParcelNames);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='ParcelFC';
%Allows you to set name for this particular analysis
if isempty(AnalysisName)
    AnalysisName=uiEnterName(['ParcelFC_',strrep(ParcelTCName,'ParcelTC_','')],['Enter name for ',AnalysisType,newline,'analysis below:']);
end

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end

if Overwrite==1 || ~isfield(AnalysisParams,'RunDate')
    AnalysisParams.RunDate=genDateString;
end
AnalysisParams.AnalysisType=AnalysisType;
AnalysisParams.AnalysisName=AnalysisName;
AnalysisParams.Parcellation=ParcelNames;
AnalysisParams.ParcelTCName=ParcelTCName;

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
    descript1='desc-parcelFCs'; %set file description name
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
            LoadPath_ParcelTC=filePaths_ParcelTC.(ParcelName){dataInd,1};
        catch
            LoadPath_ParcelTC=[];
        end
        %Pull base timecourse data. If it doesn't exist, skip run.
        if isempty(LoadPath_ParcelTC)
            continue
        end
        try
            TempLoadData = load(LoadPath_ParcelTC,'ParcelTC','AnalysisParameters');
            ParcelTC=table2array(TempLoadData.ParcelTC); 
            AnalysisParams.ParcelTCParams=TempLoadData.AnalysisParameters;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_ParcelTC]);
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
        [FCMat,FCVert] = ComputeFC(ParcelTC);  
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
        
        save(SaveNames{1,1},'FC','FCMat','FCMat_fZ','FCMat_ZNorm','FCMat_fZ_ZNorm','fZ_Degree','fZ_WholeCM','fZ_ZNorm_Degree','fZ_ZNorm_WholeCM');
    end     
toc    
end
for parcelNum=1:numParcels
    ParcelName=ParcelNames{parcelNum,1};
    parcelName=ParcelName;
  
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
    [FCMats] = CompileND(ExperimentsDir,fmriprep_table,...
        'AnalysisType',AnalysisType,...
        'LoadVarName','FCMat_fZ',...
        'DataType','ByParcellation',...
        'LoadVarFormat','Table',...
        'AnalysisName',AnalysisName,...
        'TableOrArray','Table',...
        'SubjectOrRun',SubjectOrRun,...
        'ParcelName',ParcelName);
    [Group_FC,Group_FCMat_Ts,Group_FCMat_Means,FC_Degree_Table,ReliabilityResults] = GroupMatParcellationSummaryFigures(FCMats,GroupBrainMapsDir,GroupFiguresDir,ParcelName,UseMask,UseLabels,...
        'matName','FCz','MakeFigs',MakeFigs,'ComputeReliability',ComputeReliability,'MakeBrainMaps',MakeBrainMaps);
    if ComputeReliability == 1
        FC_SplitHalfReliability=ReliabilityResults.Mat_SplitHalfReliability;
        FC_losoReliability=ReliabilityResults.Mat_losoReliability;
        save([GroupAnalysisDir,'Group_ParcelReliability'],'FC_SplitHalfReliability','FC_losoReliability');     
    end    

    save([GroupAnalysisDir,'AnalysisParams_FC'],'AnalysisParams');
    save([GroupAnalysisDir,'Group_FC'],'Group_FC');
    save([GroupAnalysisDir,'Group_FCMat'],'Group_FCMat_Ts','Group_FCMat_Means','FC_Degree_Table');  
      
end

end