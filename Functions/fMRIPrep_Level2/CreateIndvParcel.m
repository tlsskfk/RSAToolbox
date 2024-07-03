function [fmriprep_table] = CreateIndvParcel(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Template function for data processing from the fmriprep_table

%% Set default values using variable setter function
%Overwrite previously saved files (default is no or 0; yes = 1)
[Overwrite] = VariableSetter('Overwrite',0,varargin);
%Set analysis name. Default will prompt a request.
[AnalysisName] = VariableSetter('AnalysisName',[],varargin);
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Run',varargin);
BrainSpace = VariableSetter( 'BrainSpace',[],varargin);
VoxelSize = VariableSetter( 'VoxelSize',[],varargin);
ItemGroups = VariableSetter( 'ItemGroups',[],varargin);
NewVarOrder = VariableSetter( 'NewVarOrder',[],varargin);
Labels = VariableSetter( 'Labels',[],varargin);
CatParcels = VariableSetter( 'CatParcels',[],varargin);
ParcelSaveName = VariableSetter( 'ParcelSaveName',[],varargin);
[fmriprep_table_name] = VariableSetter('fmriprep_table_name',[],varargin);
MaskParams=VariableSetter( 'MaskParams',[],varargin);

%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);
if isempty(fmriprep_table_name)
[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
end
%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.
AnalysisType='IndvParcels';
%Sets analysis name as boldTC_(VoxelSize)

if ~isempty(MaskParams)
    VoxelSize=MaskParams.VoxelSize;
    ItemGroups=MaskParams.ItemGroups;
    NewVarOrder=MaskParams.NewVarOrder;
    CatParcels=MaskParams.CatParcels;
    Labels=MaskParams.Labels;
    ParcelSaveName=MaskParams.ParcelSaveName;
end
    
if isempty(VoxelSize)
    SingleSelect=1;
    VoxelSize=uiNameSelect({'3mm','2mm','1mm'},'Select target voxel size (in mm^3).',SingleSelect);
end
if isempty(BrainSpace)
    SingleSelect=1;
    BrainSpace=uiNameSelect({'MNI','TAL','Orig'},'Select brain normalization space.',SingleSelect);
end
if isempty(AnalysisName)
    AnalysisName=['UseParcels_',VoxelSize,BrainSpace];
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
if isempty(CatParcels)
    SingleSelect=0;
    CatParcels=uiNameSelect({'Yes','No','Auto'},'Combine ROIs in parcellation?',SingleSelect);
end
[filePaths,~,fsAnalysisName]=BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','AnalysisType',AnalysisType);
IndvParcelNames=filePaths.Properties.VariableNames;
UseIndvParcels=uiNameSelect(IndvParcelNames,'Select sources for Individual ROIs');
if ~iscell(UseIndvParcels)
    UseIndvParcels={UseIndvParcels};
end
AllLabels=cell(1,length(UseIndvParcels));
if contains(fsAnalysisName,'FSParcels')
    load('Parcellations/fsParcelLabels.mat','fsParcelLabels');
    for i = 1:length(UseIndvParcels)
        tempNames=uiNameSelect(fsParcelLabels.(UseIndvParcels{i})(:,2),['Select ROIs from ',UseIndvParcels{i}]);
        if ~iscell(tempNames)
            tempNames={tempNames};
        end        
        AllLabels{1,i}=tempNames;
    end
else    
    for i = 1:length(UseIndvParcels)
        UseParcelName=UseIndvParcels{i};
        for f = 1:height(filePaths)
            try
                tempload=load(filePaths.(UseParcelName){f,1},'UseLabels');
                tempNames=uiNameSelect(tempload.UseLabels,['Select ROIs from ',UseIndvParcels{i}]);
                break
            catch
                continue
            end
        end
        if ~iscell(tempNames)
            tempNames={tempNames};
        end        
        AllLabels{1,i}=tempNames;
    end
end
%
[~,AllGroupParcelNames]=getFolderAndFileNames('Parcellations/',[VoxelSize,BrainSpace]);
GroupParcelNames=uiNameSelect([{'None'};AllGroupParcelNames(:)],'Select group parcels to use');  
if ~strcmpi(GroupParcelNames{1,1},'None')
GroupLabels=cell(1,length(GroupParcelNames));
for i = 1:length(GroupParcelNames)
    load(['Parcellations/',GroupParcelNames{i}],'UseLabels');
    tempNames=uiNameSelect(UseLabels,['Select ROIs from ',GroupParcelNames{i}]);
    if ~iscell(tempNames)
        tempNames={tempNames};
    end
    GroupLabels{1,i}=tempNames;
end
else
    GroupLabels=[];
end

GroupMask=[];
GroupMaskLabels=[];
pCount=1;
if isempty(ParcelSaveName)
    ParcelSaveName=uiEnterName(['_',VoxelSize,BrainSpace],'Enter parcel save name');
end
MaskParams.ParcelSaveName=ParcelSaveName;
if ~isempty(GroupLabels)
    for i = 1:length(GroupParcelNames)
        load(['Parcellations/',GroupParcelNames{i}],'UseLabels','UseMask');
        SelectedLabels=GroupLabels{1,i};
        if isempty(GroupMask)
            GroupMask=UseMask*0;
        end
        for j = 1:length(SelectedLabels)
            ROIName=SelectedLabels{j};
            ROINum=find(ismember(UseLabels,ROIName));
            GroupMask(UseMask==ROINum)=pCount;
            if any(ismember(GroupMaskLabels,ROIName))
                stopSwitch=0;
                stopCt=1;
                while stopSwitch == 0
                    if stopCt == 1
                        ROIName=[ROIName,'Grp1'];
                    else
                        ROIName(1,end)=num2str(stopCt);
                    end
                    stopCt=stopCt+1;
                    if ~any(ismember(GroupMaskLabels,ROIName))
                        stopSwitch=1;
                    end
                end
            end
            GroupMaskLabels{pCount,1}=ROIName;
            pCount=pCount+1;
        end
    end
else
    load(['Parcellations/',AllGroupParcelNames{1}],'UseMask');
    GroupMask=UseMask*0;
    GroupMaskLabels=cell(1);
end
%% BIDsTable loop: Iterate through fmriprep_table and perform analysis
AllUseLabels=GroupMaskLabels;

MakeLabels=0;

AllMasks=cell(height(fmriprep_table),1);
AllLabelsIni=[];
for sNum = 1:height(filePaths)
    spCount=pCount;
    IndvMask=GroupMask;
    
    try
        load(filePaths.(UseIndvParcels{1}){sNum,1},'UseMask','UseLabels');
    catch   
        continue
    end
    for i = 1:length(UseIndvParcels)
        IndvPName=UseIndvParcels{i};
        SelectedLabels=AllLabels{1,i};
        load(filePaths.(IndvPName){sNum,1},'UseMask','UseLabels');
        for j = 1:length(SelectedLabels)
            ROIName=SelectedLabels{j};
            ROINum=find(ismember(UseLabels,ROIName));
            IndvMask(UseMask==ROINum)=spCount;
            if~isempty(AllUseLabels{1,1})
                if MakeLabels==0
                    if any(ismember(AllUseLabels,ROIName))
                        stopSwitch=0;
                        stopCt=1;
                        while stopSwitch == 0
                            if stopCt == 1
                                ROIName=[ROIName,'Indv1'];
                            else
                                ROIName(1,end)=num2str(stopCt);
                            end
                            stopCt=stopCt+1;
                            if ~any(ismember(AllUseLabels,ROIName))
                                stopSwitch=1;
                            end
                        end
                    end
                    AllUseLabels{spCount,1}=ROIName;
                end   
            elseif MakeLabels==0
                AllUseLabels{spCount,1}=ROIName;
            end
            spCount=spCount+1;
        end      
    end
    if isempty(AllLabelsIni)
        AllLabelsIni=AllUseLabels;
        AllLabelsfin=AllUseLabels;
    end
    if any(ismember(CatParcels,'Yes'))
        if isempty(ItemGroups)
            [~,Labels,ItemGroups] = uiParcelCat(AllUseLabels,'Combine ROIs for parcel');                   
        end
        AllUseLabels=Labels(:);
        IndvMask2=IndvMask*0;
        for i = 1:length(ItemGroups)
            ItemGroup=ItemGroups{1,i};
            for j = 1:length(ItemGroup)
                ROINum=find(ismember(AllLabelsIni,ItemGroup{1,j}));
                IndvMask2(IndvMask==ROINum)=i;
            end
        end

        IndvMask=IndvMask2;       
    elseif any(ismember(CatParcels,'Auto'))
        IndvMask2=IndvMask*0;
        AllLabelsfin=strrepCell(AllLabelsIni,'L_','');
        AllLabelsfin=strrepCell(AllLabelsfin,'R_','');
        AllUseLabels=unique(AllLabelsfin);
        for m= 1:length(AllUseLabels)
            useNums=find(ismember(AllLabelsfin,AllUseLabels{m}));
            IndvMask2(ismember(IndvMask,useNums))=m;
        end
        IndvMask=IndvMask2;
    end          
    if MakeLabels==0 
        if isempty(NewVarOrder)
            [NewVarOrder] = uiNameSelect(AllUseLabels);
        end
        [~,NewVarNum]=ismember(AllUseLabels,NewVarOrder); 
        AllUseLabels=NewVarOrder;
        MakeLabels=1;
    end
    IndvMask2=IndvMask*0;
    for k = 1:length(NewVarNum)
        IndvMask2(IndvMask==k)=NewVarNum(k,1);
    end
    AllMasks{sNum,1}=IndvMask2;
end
UseLabels=AllUseLabels;
iniPercentComplete=0;
for dataInd=useIndicies    
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. ',AnalysisType,' - ',AnalysisName]);
        iniPercentComplete=PercentComplete;
    end
    
    %% Set save directory and save name
    SaveDir=strrep(fmriprep_table.funcDir{dataInd,1},'/func/',['/',AnalysisType,'/',AnalysisName,'/',ParcelSaveName,'/']);
    SaveDir=strrep(SaveDir,'/fmriprep/','/matlab/');
    SaveName = fmriprep_table.preproc_bold{dataInd,1};
    SaveName=strrep(SaveName,'.nii','');
    SaveName=strrep(SaveName,'.gz','');
    if fmriprep_table.numRuns(dataInd,1)>1 && bySS==1
        SaveName=strrep(SaveName,'_run-01','');
    end
    descript1='desc-IndvParcel'; %set file description name
    SaveName=strrep(SaveName,'desc-preproc_bold',descript1);
    
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    SavePrefix=[ExperimentsDir,SaveDir];    
    if ~exist(SavePrefix,'file')
        mkdir(SavePrefix);
    end    
    SaveNames{1,1}=[SavePrefix,SaveName,'.mat'];
    if exist(SaveNames{1,1},'file')~=0 && Overwrite==0
        disp(['Skipping-- file exists: ',SaveName]);
        continue
    end          
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;    
        UseMask=AllMasks{loadInd,1};
        save(SaveNames{1,1},'UseMask','UseLabels');     
    end     
   toc
end 
GroupSaveDir=['Parcellations/IndvParcels/',fmriprep_table_name,'/'];
if ~exist(GroupSaveDir,'file')
    mkdir(GroupSaveDir);
end    
MaskParams.VoxelSize=VoxelSize;
MaskParams.ItemGroups=ItemGroups;
MaskParams.NewVarOrder=NewVarOrder;
MaskParams.CatParcels=CatParcels;
MaskParams.Labels=Labels;
MaskParams.BrainSpace=BrainSpace;
MaskParams.ParcelSaveName=ParcelSaveName;

MasksBySS=AllMasks(fmriprep_table.numRuns_bySub~=0,:);
StatsTable=fmriprep_table(fmriprep_table.numRuns_bySub~=0,{'sub'});
emptySS=cellfun(@isempty, MasksBySS);
MissingSS=StatsTable.sub(emptySS,:);
MasksBySS(emptySS,:)=[];
StatsTable(emptySS,:)=[];

ROIStats=nan(height(StatsTable),length(UseLabels));
for i = 1:height(StatsTable)
    tempInd=unique(MasksBySS{i,1}(:));
    tempInd(1,:)=[];
    tempCts=histcounts(MasksBySS{i,1}(:),'BinMethod','integers');
    try
    ROIStats(i,tempInd)=tempCts(1,2:end);
    end
end
StatsTable=[StatsTable,array2table(ROIStats,'VariableNames',UseLabels)];

MasksMat=cell2nDMAT(MasksBySS);
MaskByParcel=repmat(UseMask*0,[1,1,1,length(UseLabels)]);
for i = 1:length(UseLabels)
    MaskByParcel(:,:,:,i) = sum(single(MasksMat==i),4);
end
MaskByParcel(isnan(MaskByParcel))=0;
temp=cat(4,MaskByParcel(:,:,:,1)*0,MaskByParcel);
[~,UseMask]=max(temp,[],4);
UseMask=UseMask-1;
save([GroupSaveDir,ParcelSaveName],'MaskParams','StatsTable','MaskByParcel','UseMask','UseLabels','MissingSS');


