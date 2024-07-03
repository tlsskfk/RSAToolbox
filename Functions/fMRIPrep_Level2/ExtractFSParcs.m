function [fmriprep_table] = ExtractFSParcs(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Template function for data processing from the fmriprep_table

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Run',varargin);
fsDir=VariableSetter('fsDir',[],varargin);
fpDir=VariableSetter('fpDir',[],varargin);
WorkDir=VariableSetter('WorkDir',[],varargin);
Resample = VariableSetter( 'Resample',[],varargin);
VoxelSize = VariableSetter( 'VoxelSize',[],varargin);


%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='IndvParcels';
%Sets analysis name as boldTC_(VoxelSize)
if isempty(VoxelSize)
    SingleSelect=1;
    VoxelSize=uiNameSelect({'3mm','2mm','1mm'},'Select target voxel size (in mm^3).',SingleSelect);
end
AnalysisName=['FSParcels',VoxelSize];
if isempty(fsDir)
    fsDir=uigetdir(ExperimentsDir,'Select FreeSurfer directory');
    fsDir=[fsDir,'/'];
    fsDir=strrep(fsDir,'\','/');
end
if isempty(fpDir)
    fpDir=uigetdir(ExperimentsDir,'Select fmriPrep directory');
    fpDir=[fpDir,'/'];
    fpDir=strrep(fpDir,'\','/');
end
if isempty(WorkDir)
    WorkDir=uigetdir('/Users/','Select work directory');
    WorkDir=[WorkDir,'/'];
    WorkDir=strrep(WorkDir,'\','/');
end
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

%Set Error (or other variable) table column in case variable doesnt exist.
%Error column is used for recording errors picked up using try, catch,
%operators
if sum(ismember(fmriprep_table.Properties.VariableNames,'Error'))==0
    fmriprep_table.Error=cell(TotalRuns,1);
end
DirLetter=lower(WorkDir(1,1));
wslWorkDir=strrep(WorkDir,WorkDir(1,1:2),['/mnt/',DirLetter]);  

%% BIDsTable loop: Iterate through fmriprep_table and perform analysis
iniPercentComplete=0; %Used to display progress
LoadVars=[];
load('Parcellations/fsParcelLabels.mat','fsParcelLabels');
sfParcName{1,1}='aparc.DKTatlas+aseg.mgz';
sfParcName{2,1}='aparc.a2009s+aseg.mgz';
sfParcName{3,1}='aparc+aseg.mgz';
sfParcName{4,1}='aseg.auto_noCCseg.mgz';
talFileName='transforms/talairach.xfm';
for dataInd=useIndicies    
    tic
    delete([WorkDir,'*']);
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    
    %% Set save directory and save name
    SaveDir=strrep(fmriprep_table.funcDir{dataInd,1},'/func/',['/',AnalysisType,'/',AnalysisName,'/']);
    SaveDir=strrep(SaveDir,'/fmriprep/','/matlab/');
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-fsParcels'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        try
        loadInd=dataInd+run-1;
        SSfsDir=[fsDir,'sub-',fmriprep_table.sub{loadInd,1},'/mri/'];
        %skip previous errors
%         if ~isempty(fmriprep_table.Error{loadInd,1})
%             continue
%         end    
        %% set load paths and variable names
        brainMaskName=fmriprep_table.GM_probseg{dataInd,1};
        copyfile([fpDir,'sub-',fmriprep_table.sub{loadInd,1},'/anat/',brainMaskName],WorkDir,'f');
        copyfile([SSfsDir,talFileName],WorkDir,'f');
        for fNum=1:length(sfParcName)         
            copyfile([SSfsDir,sfParcName{fNum,1}],WorkDir,'f');           
        end
        
        for r = 1:4
            baseName=strrep(sfParcName{r,1},'.mgz','');
            runCode = [];
            runCode=[runCode,'cd',newline];
            runCode=[runCode,'export PATH=/usr/local/go/bin:$PATH',newline];
            runCode=[runCode,'export FS_LICENSE=/home/david/license.txt',newline];
            runCode=[runCode,'export FREESURFER_HOME=/usr/local/freesurfer/7.3.2',newline];
            runCode=[runCode,'export FSFAST_HOME=/usr/local/freesurfer/7.3.2/fsfast',newline];
            runCode=[runCode,'source /usr/local/freesurfer/7.3.2/SetUpFreeSurfer.sh',newline];
            runCode=[runCode,'export SUBJECTS_DIR=/.../data/ftadel/freesurfer',newline];
            runCode=[runCode,'export FUNCTIONALS_DIR=/.../data/ftadel/freesurfer',newline];
            runCode=[runCode,'source $FREESURFER_HOME/FreeSurferEnv.sh',newline];
            runCode=[runCode,'cd ',wslWorkDir,newline];
            runCode=[runCode,'mri_convert ',baseName,'.mgz ',baseName,'.tlrc.mgz -at talairach.xfm -rt nearest -oc 0 0 0',newline];
            runCode=[runCode,'mri_label2vol --reg $FREESURFER_HOME/average/mni152.register.dat --seg ',baseName,'.tlrc.mgz  --temp ',brainMaskName,' --o ',baseName,'.mni152.nii',newline];
            fsCodeFileName=[WorkDir,'fsCode.sh'];
            WSLfsCodeFileName=[wslWorkDir,'fsCode.sh'];
            fileID = fopen(fsCodeFileName,'w');
            fprintf(fileID,'%s',runCode);
            fclose(fileID);
            [~,temp]=system(['wsl bash -c ',WSLfsCodeFileName]);
            brainMap=load_nii([WorkDir,baseName,'.mni152.nii']);
            baseName=strrep(baseName,'.','');
            baseName=strrep(baseName,'+','');         
            tempParcel=single(brainMap.img);
            if ~isempty(Resample)
                tempParcel=imresize3(tempParcel,Resample,'nearest');
            end
            
            [UseMask,UseLabels]=Parcel2Mask(tempParcel,fsParcelLabels.(baseName));
            SavePrefix=[ExperimentsDir,SaveDir,baseName,'/'];    
            if ~exist(SavePrefix,'file')
                mkdir(SavePrefix);
            end    
            SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
            if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
                disp(['Skipping-- file exists: ',SaveName]);
                continue
            end             
            save(SaveNames{1,1},'UseMask','UseLabels'); 
        end
            catch
        disp('Error')
        end    
    end     
    delete([WorkDir,'*']);
    toc
end
end 

function [UseMask,UseLabels]=Parcel2Mask(Parcel,Key)
UseMask=Parcel*0;
UseLabels=Key(:,2);
for i = 1:size(Key,1)
    UseMask(Parcel==Key{i,1})=i;
end
end