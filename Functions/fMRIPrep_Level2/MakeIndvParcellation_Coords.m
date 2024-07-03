function [AnalysisParameters] = MakeIndvParcellation_Coords(fmriprep_table,ExperimentsDir,varargin)
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
[AnalysisName] = VariableSetter('AnalysisName',[''],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun',[],varargin);
[ROInames] = VariableSetter('ROInames',[],varargin);
[ROIshape] = VariableSetter('ROIshape',[],varargin);
[ROIradius] = VariableSetter('ROIradius',[],varargin);
[loadCoords] = VariableSetter('loadCoords',[],varargin);
% Apply graymatter probability mask. Enter an number between 0 (most
% lenient) to 1 (most conservative). Default is empty (no mask)
[gmMask] = VariableSetter('gmMask',[],varargin);
%Input predefined analysis parameters analysis parameters
[AnalysisParameters] = VariableSetter('AnalysisParameters',[],varargin);

ExperimentsDir=strrep(ExperimentsDir,'\','/');
GroupDir=strrep(ExperimentsDir,'Experiments/','GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
FigureDir=strrep(GroupDir,'/GroupAnalysis/','/GroupFigures/');
BrainMapDir=strrep(GroupDir,'/GroupAnalysis/','/GroupBrainMaps/');

MaxSphereSLSize=[1,7;2,33;3,123;4,257;5,515;6,925;7,1419;8,2109;9,3017;10,4169;11,5575;12,7153;13,9171;14,11513;15,14147];
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
AnalysisType='CoordParcels';
%Allows you to set name for this particular analysis

if isempty(fmriprep_table_name)
    [~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
%% Compile filepaths for input files for the analysis
[filePaths_BaseTC,~,BaseTCName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','func','AnalysisName','raw3mm','TitleTextName','Select fMRI input for GLM:');
filePaths_BaseTC=filePaths_BaseTC.(BaseTCName);
brainDims=[];
for i = 1:length(filePaths_BaseTC)
    try
        tempLoad=load(filePaths_BaseTC{i,1},'brain_mask');
        brainDims=size(tempLoad.brain_mask);
    catch
        continue
    end
    break
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

if isempty(loadCoords)
    SingleSelect=1; %Allows only a single value to be selected.
    [loadCoords] = uiNameSelect({'Yes','No'},'Load saved ROI coordinates?',SingleSelect);
    if strcmpi(loadCoords,'Yes')
        loadCoords=1;
    else
        loadCoords=0;
    end
end    
if loadCoords==1
    [filePaths_loadCoords,~,loadCoordNames] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Subject','AnalysisType','ROICoords','TitleTextName','Select ROI coordinates to load:');
    filePaths_loadCoords=filePaths_loadCoords.(loadCoordNames); 
    loadROINames=[];
    for i = useIndicies
        if ~isempty(filePaths_loadCoords{i,1})
            load(filePaths_loadCoords{i,1},'ROICoords');
            loadROINames=unique([loadROINames;ROICoords.Properties.VariableNames']);
        end
    end
    loadCoordSuffix=['_',loadCoordNames];
else
    filePaths_loadCoords=[];
    loadCoordSuffix='';
end
%% Apply graymatter mask
gmAffix=[];
if isempty(gmMask)
    gmMask=uiEnterName('',['Apply graymatter mask? Leave blank for no.',newline,'Enter value between 0 and 1', newline,'0 = most lenient; 1 is most conservative']);    
    if ~isempty(gmMask)
        gmMask=str2num(gmMask);
        gmAffix=['_gm',num2str4filename(gmMask,2)];
    end
end

if ~isempty(gmMask)
    [filePaths_gmMask,~,gmMaskName] = BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType','anat','TitleTextName','Select anat data for graymatter mask:');
    filePaths_gmMask=filePaths_gmMask.(gmMaskName);
end

%% Set searchligh setting
if isempty(ROIshape)
    SingleSelect=1; %Allows only a single value to be selected.
    [ROIshape] = uiNameSelect({'Sphere','Cube'},'Select ROI shape:',SingleSelect);
end

if isempty(ROIradius)
    ROIradiusName=uiEnterName('3','Enter ROI radius');
    ROIradius=str2num(ROIradiusName);
else
    ROIradiusName=num2str(ROIradius);
end

if isempty(ROInames)
    numROIsName=uiEnterName('3','Enter number of ROIs');
    numROIs=str2num(numROIsName);
    ROInames=cell(numROIs,2);
    for i = 1:numROIs
        ROInames{i,1}=uiEnterName(['ROI',num2str(i)],['Enter name for ROI # ',num2str(i)]);
        if loadCoords == 1
            SingleSelect=0; 
            [ROInames{i,2}] = uiNameSelect(loadROINames,'Select coords to include:',SingleSelect);
            loadROINames(ismember(loadROINames,ROInames{i,2}),:)=[];
        end
    end
elseif isempty(ROInames{1,2})
    numROIs=str2num(numROIsName);
    ROInames=cell(numROIs,2);
    for i = 1:numROIs
        if loadCoords == 1
            SingleSelect=0; 
            [ROInames{i,2}] = uiNameSelect(loadROINames,'Select coords for ',ROInames{i,1},' to include:',SingleSelect);
            loadROINames(ismember(loadROINames,ROInames{i,2}),:)=[];
        end
    end

else    
   numROIs=size(ROInames,1);
   numROIsName=num2str(numROIs);
end

if isempty(AnalysisName)
    AnalysisName=uiEnterName(['CoordParcel_Shape-',ROIshape,'_Rad-',ROIradiusName,'_',numROIsName,'ROIs',gmAffix,loadCoordSuffix],['Enter name for ',AnalysisType,newline,'analysis below:']);
    if length(AnalysisName) > namelengthmax
        AnalysisName=uiEnterName(AnalysisName(1,1:namelengthmax),['Please shorten name: max length = ',num2str(namelengthmax)]);
    end
end
AnalysisParameters.ROInames=ROInames;
AnalysisParameters.ROIradius=ROIradius;
AnalysisParameters.ROIshape=ROIshape;
AnalysisParameters.SubjectOrRun=SubjectOrRun;
AnalysisParameters.gmMask=gmMask;

iniPercentComplete=0; %Used to display progress
AllMasks=[];
if bySS == 1
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'_bySS/',fmriprep_table_name,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'_bySS/',fmriprep_table_name,'/'];
else
    GroupAnalysisDir=[GroupDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',fmriprep_table_name,'/'];
    GroupBrainMapsDir=[BrainMapDir,fmriprep_table_name,'/',AnalysisType,'/',AnalysisName,'/',fmriprep_table_name,'/'];
end
if ~exist(GroupAnalysisDir,'file')
    mkdir(GroupAnalysisDir);     
end
if ~exist(GroupBrainMapsDir,'file')
    mkdir(GroupBrainMapsDir);     
end   
ROIInfo=[];
AllMasksByROI=zeros(brainDims(1,1),brainDims(1,2),brainDims(1,3),size(ROInames,1),length(useIndicies));
AllMasksBySS=zeros(brainDims(1,1),brainDims(1,2),brainDims(1,3),length(useIndicies));
counter=1;
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
        SaveName=strrep(SaveName,'_run-01','');
        SaveName=strrep(SaveName,'_run-1','');
    end
    descript1='desc-CoordParcel'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    %% If analysis involves parcellations:
    SaveNames=cell(1);
    SavePrefix=[ExperimentsDir,SaveDir,'/']; 
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    else  
        SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
        if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
            disp(['Skipping-- all files exist: ',SaveName]);
        end
    end
    %% Determine number of runs
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
        rowID=fmriprep_table.sub{dataInd,1};
    else
        numRuns=1;
        rowID=[fmriprep_table.sub{dataInd,1},'_run',num2str(fmriprep_table.run(dataInd,1))];
    end
    count=1;
    %initialize variable names
    BaseTC=cell(1);
    brain_mask=cell(1);
    Use_gmMask=[];    
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;
        %skip previous errors
        if ~isempty(fmriprep_table.Error{loadInd,1})
            continue
        end    
        %% set load paths and variable names
        %pull load paths
        try
            LoadPath_BaseTC=filePaths_BaseTC{loadInd,1};
        catch
            LoadPath_BaseTC=[];
        end     
        %Pull base timecourse data. If it doesn't exist, skip run.
        try
            variableInfo = who('-file', LoadPath_BaseTC);
            if any(ismember(variableInfo,'brain_mask'))                
                TempLoadData = load(LoadPath_BaseTC,'brain_mask'); 
            end
            brain_mask{count,1}=TempLoadData.brain_mask;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadPath_BaseTC]);
            continue
        end
        if ~isempty(gmMask) && isempty(Use_gmMask)
%             try
%                 TempLoadData = load(filePaths_gmMask{loadInd,1},'WM_probseg');
%                 Use_wmMask=single(TempLoadData.WM_probseg < 0.99);
%             catch
%                 disp(['White matter mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
%                 Use_wmMask=(brain_mask{count,1}*0)+1;
%             end
% 
%             try
%                 TempLoadData = load(filePaths_gmMask{loadInd,1},'CSF_probseg');
%                 Use_csfMask=single(TempLoadData.CSF_probseg < 0.25);
%             catch
%                 disp(['CSF mask error-- mask not applied',filePaths_gmMask{loadInd,1}]);
%                 Use_csfMask=(brain_mask{count,1}*0)+1;
%             end            
            try
                LoadPath_gmMask=filePaths_gmMask{loadInd,1}; %GM_probseg
                TempLoadData = load(LoadPath_gmMask,'GM_probseg');
                Use_gmMask=single(TempLoadData.GM_probseg > gmMask);
            catch
                disp(['Gray matter mask error-- mask not applied',LoadPath_gmMask]);
                Use_gmMask=(brain_mask{count,1}*0)+1;
            end
        elseif isempty(Use_gmMask)
            Use_gmMask=(brain_mask{count,1}*0)+1;
%             Use_wmMask=(brain_mask{count,1}*0)+1;
%             Use_csfMask=(brain_mask{count,1}*0)+1;
        end          
        count=count+1;
    end     
    if count==1
        disp(['Skipping subject-- no input files or variables exist']);
        continue
    end  
    brain_mask=single(sum(single(cell2nDMAT(brain_mask)~=0),4)==size(brain_mask,1)).*Use_gmMask;
    if loadCoords==0
        [UseLabels,Coords] = uiROICoordAssigment(ROInames(:,1),['Enter coords for ',SaveName]);
        coordFilter=single(sum(abs(Coords),2))==0;
        Coords=num2cell(Coords,2);
        ROICoords=cell2table(Coords','VariableNames',UseLabels);
    else
        load(filePaths_loadCoords{dataInd,1},'ROICoords');
        UseLabels=ROInames(:,1);
        Coords=cell(length(UseLabels),1);
        coordFilter=zeros(length(UseLabels),1);
        for i = 1:length(UseLabels)
            try
                Coords{i,1}=cell2mat(table2cell(ROICoords(:,ROInames{i,2}))');
            catch
                disp('No ROI')
                coordFilter(i,1)=1;
            end
        end
    end
    UseMask=brain_mask*0;
    ROIsizes=zeros(1,length(UseLabels));
    for roiNum=1:length(UseLabels)
        if coordFilter(roiNum,1)==0
            [MatCoords] = CoordConvert_MNI2Mat(Coords{roiNum,1},size(brain_mask));
            spherecoords=[];
            for i = 1:size(MatCoords,1)
                [ tempCoords ] = coord2spherecoords( MatCoords(i,:),ROIradius );
                spherecoords=[spherecoords;tempCoords];
            end
            [spherecoords]=filtercoords(spherecoords,size(UseMask,1),size(UseMask,2),size(UseMask,3),1);
            tempMat=coords2mat(spherecoords,UseMask,ones(size(spherecoords,1),1));
            tempMat=tempMat.*brain_mask;
            ROIsizes(1,roiNum)=sum(tempMat(:));
            AllMasksByROI(:,:,:,roiNum,counter)=tempMat;
            UseMask(tempMat==1)=roiNum;
        end
    end
    AllMasksBySS(:,:,:,counter)=UseMask;
    if isempty(ROIInfo)
        ROIInfo=array2table(ROIsizes,'VariableNames',UseLabels,'RowNames',{rowID});
    else
        ROIInfo=[ROIInfo;array2table(ROIsizes,'VariableNames',UseLabels,'RowNames',{rowID})];
    end
    save(SaveNames{1,1},'UseMask','UseLabels','AnalysisParameters','ROInames','ROICoords');  
    counter=counter+1;
end

GroupAllMasks=nansum(AllMasksByROI,5);
SaveBrik_3mmMNI(GroupAllMasks,ROInames(:,1),[GroupBrainMapsDir,'ROIOverlapMap']);
SaveBrik_3mmMNI(AllMasksBySS,ROIInfo.Properties.RowNames(:),[GroupBrainMapsDir,'AllSS_Map']);
save([GroupAnalysisDir,'MapInfo'],'AllMasksByROI','GroupAllMasks','ROIInfo','AnalysisParameters');
GroupSaveDir=['Parcellations/IndvCoords/',fmriprep_table_name,'/'];
if ~exist(GroupSaveDir,'file')
    mkdir(GroupSaveDir);
end  
GroupAllMasks(isnan(GroupAllMasks))=0;
temp=cat(4,GroupAllMasks(:,:,:,1)*0,GroupAllMasks);
[~,UseMask]=max(temp,[],4);
UseMask=UseMask-1;
save([GroupSaveDir,AnalysisName],'UseMask','UseLabels');

