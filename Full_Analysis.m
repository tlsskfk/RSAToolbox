clear
close all
%fmriprep_table_name='FaceDP_20042002_noDPs';
%ExperimentsDir=['J:\fMRI\LettersAndDigits\Experiments\'];
%ExperimentsDir='F:\Data\fMRI\WordReading\Experiments\';
ExperimentsDir='G:\fMRI\gradCPT\Experiments\';
ExpName='/gradCPT/';
fmriprep_table_name=['DoDgradCPT'];
try
    load(['fmriprep_table',ExpName,fmriprep_table_name]);
catch
    fmriprep_table=[];
end
gradCPT=1;
if isempty(ExperimentsDir)
    ExperimentsDir=uigetdir('/');
    ExperimentsDir=[ExperimentsDir,'/'];
end

if isempty(fmriprep_table)
    disp(['Setting up fmriprep_table for subsequent analyses']);
    fmriprep_table=Create_fMRIPrep_Table(ExperimentsDir);
    save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');
end

% %%% Updating old data to BIDs format %%%
runID=uiEnterName('0','Add gradCPT OrigIDs; 1: Yes 0:No');
if str2num(runID)==1
    disp(['Adding original ID info']);
    [fmriprep_table] = fmriprep_table_AddOrigIDs(ExperimentsDir,fmriprep_table);
    save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');
end

%Extract bold timecourses
disp(['Converting raw fmriprep BOLD output to matlab while applying minor cleaning']);
VoxelSize=[];
if isempty(VoxelSize)
    SingleSelect=1;
    VoxelSize=uiNameSelect({'3mm','2mm','1mm'},'Select target voxel size (in mm^3).',SingleSelect);
end
[fmriprep_table] = ExtractAndCleanBoldTCs(ExperimentsDir,fmriprep_table,'Resample',[61,73,61],'VoxelSize',VoxelSize,'Overwrite',0);
[fmriprep_table] = ExtractAnatMasks(ExperimentsDir,fmriprep_table,'Resample',[61,73,61],'VoxelSize',VoxelSize,'Overwrite',0);
save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');


%%%Pull out fmriprep confound TCs, compute motion stats
% disp(['Extracting fmriprep confound TCs']);
[ExperimentsDir,fmriprep_table] = ComputeConfoundTC(ExperimentsDir,fmriprep_table,'Overwrite',1);
save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');
disp(['Computing summary fmriprep motion stats']);
[fmriprep_table] = fMRIPrep_Table_AddStartInds(fmriprep_table)
[fmriprep_table] = GenerateMotionStats(ExperimentsDir,fmriprep_table);
save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');

%Extract and compute gradCPT specific variables
if gradCPT==1
    [fmriprep_table] = gradCPT_BehaviorDataMigrate(fmriprep_table,ExperimentsDir);
    disp(['Extracting gradCPT Events']);
    [ fmriprep_table ] = gradCPT_ExtractEvents(fmriprep_table,ExperimentsDir,'Overwrite',0);
    save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');   
    disp(['Computing gradCPT behavioral vars by subject']);
    [fmriprep_table] = ComputeBehavioralVars(fmriprep_table,ExperimentsDir,'Overwrite',0,'SubjectOrRun','Subject');
    disp(['Computing gradCPT behavioral vars by run']);
    [fmriprep_table] = ComputeBehavioralVars(fmriprep_table,ExperimentsDir,'Overwrite',0,'SubjectOrRun','Run');    
else
    [fmriprep_table] = LoadRawBehavior(fmriprep_table,ExperimentsDir);
    save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');
    [ExperimentsDir,fmriprep_table,fmriprep_table_name] = ComputeBehaviorTCs(fmriprep_table,ExperimentsDir);
    save(['fmriprep_table',ExpName,fmriprep_table_name],'fmriprep_table');
end
%Optional level 1 GLM-based cleaning
%[BIDsTable] = BIDsComputeResids_byRun(BIDsTable,ExpDir,varargin)
%save(['BIDsTables/',BIDsTableName],'BIDsTable','ExpDir');  
% [ fmriprep_table ] = GLM_ComputeResids(fmriprep_table,ExperimentsDir,'Overwrite',0);
% 
% [ fmriprep_table ] = GLM_ComputeTs(fmriprep_table,ExperimentsDir,'Overwrite',0);
% 
% [ fmriprep_table ] = ComputeParcellationTCs(fmriprep_table,ExperimentsDir,'Overwrite',0);

%[ fmriprep_table ] = GLM_ComputeRSMs(fmriprep_table,ExperimentsDir,'Overwrite',0);


