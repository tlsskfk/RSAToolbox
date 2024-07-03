function BrainMap2GroupParcel(varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Set resample size
[Resample] = VariableSetter('Resample',[],varargin);
[BrainSpace] = VariableSetter('BrainSpace',[],varargin);
%Select reference parcellation for name
ParcelNames=cell(1);
ParcelDirInfo=dir('Parcellations/');
pcount = 1;
for i = 1:size(ParcelDirInfo,1)
    if ParcelDirInfo(i).isdir==0
        ParcelNames{pcount,1}=strrep(ParcelDirInfo(i).name,'.mat','');
        pcount=pcount+1;
    end
end

voxelOverwrite='yes';


%% Initialize variables
[editParcel] = uiNameSelect({'yes','no'},'Add to existing parcellation?',1); 
if strcmpi(editParcel,'yes')
    [voxelOverwrite] = uiNameSelect({'yes','no'},'Overwrite existing ROIs',1); 
    LoadParcelName=uiNameSelect(ParcelNames,'Select existing parcellation to edit.',1);
    tempLoadParcel=load(['Parcellations/',LoadParcelName]);
    UseMask = tempLoadParcel.UseMask;
    UseLabels = tempLoadParcel.UseLabels(:);  
    Resample=size(UseMask);
    count=size(UseLabels,1)+1;    
    if isfield(tempLoadParcel,'ParcelInfo')
        ParcelInfo=tempLoadParcel.ParcelInfo;
        if isfield(ParcelInfo,'CurrentPass')
            CurrentPass=ParcelInfo.CurrentPass+1;
        else
            CurrentPass=1;
        end
        ParcelInfo.CurrentPass=CurrentPass;
        PassStr=['Pass',num2str(CurrentPass)];
        ParcelInfo.(PassStr).Info.voxelOverwrite=voxelOverwrite;
        ParcelInfo.(PassStr).Info.LoadParcelName=LoadParcelName;
        ParcelInfo.(PassStr).Info.Resample=Resample;
    else    
        %%%%% write code here to get ParcelInfo for loaded parcel if it doesnt
        %%%%% already exist.
        %%%%%
        %%%%%
        ParcelInfo=struct;
        CurrentPass=1;
        ParcelInfo.CurrentPass=CurrentPass;
        PassStr=['Pass',num2str(CurrentPass)];
        ParcelInfo.(PassStr).Info.voxelOverwrite=voxelOverwrite;
        ParcelInfo.(PassStr).Info.LoadParcelName=LoadParcelName;
        ParcelInfo.(PassStr).Info.Resample=Resample;        
    end
else
    ParcelInfo=struct;
    CurrentPass=1;
    ParcelInfo.CurrentPass=CurrentPass;
    PassStr=['Pass',num2str(CurrentPass)];
    ParcelInfo.(PassStr).Info.LoadParcelName='None';
end
   
if isempty(BrainSpace)    
    [BrainSpace] = uiNameSelect({'mni','tal'},'Select brain space:',1); 
end
ParcelInfo.(PassStr).Info.BrainSpace=BrainSpace;
[flipMaps] = uiNameSelect({'yes','no'},'LF flip maps?',1);   
ParcelInfo.(PassStr).Info.flipMaps=flipMaps;
[fileN,dirN]=uigetfile('*.*','Select initial map (.nii .brik or .mat) to convert to parcellation');
if contains(fileN,'.nii','IgnoreCase',true)
    brainMap=load_nii([dirN,fileN]);
    brainMap=brainMap.img;
elseif contains(fileN,'.brik','IgnoreCase',true)
    [~,brainMap,info] = BrikLoad ([dirN,fileN]);
elseif contains(fileN,'.mat','IgnoreCase',true)
    SingleSelect=1; %Allows only a single value to be selected.
end
if strcmpi(flipMaps,'yes')
    brainMap=flip(brainMap,1);
end
if isempty(Resample)
    Resample = [size(brainMap,1),size(brainMap,2),size(brainMap,3)];
else    
    brainMap=imresize3(brainMap,Resample);
end
ParcelInfo.(PassStr).InitialMap=cell(1,3);
ParcelInfo.(PassStr).InitialMap{1,1}=fileN;
ParcelInfo.(PassStr).InitialMap{1,2}=brainMap;
OrigBrainMap=brainMap;
numSubtractMaps=uiEnterName('0',['Enter # subtraction-maps',newline,'to load']);
numSubtractMaps=str2num(numSubtractMaps);
if numSubtractMaps~=0
    SubtractMaps=cell(numSubtractMaps,2);
    ParcelInfo.(PassStr).SubtractMaps=cell(numSubtractMaps,3);
    for i = 1:numSubtractMaps
        [fileN,dirN]=uigetfile([dirN,'*.*'],'Select subtraction map (.nii .brik or .mat)');
        if contains(fileN,'.nii','IgnoreCase',true)
            subMap=load_nii([dirN,fileN]);
            SubtractMaps{i,1}=subMap.img;
        elseif contains(fileN,'.brik','IgnoreCase',true)
            [~,SubtractMaps{i,1},info] = BrikLoad ([dirN,fileN]);
        elseif contains(fileN,'.mat','IgnoreCase',true)
            SingleSelect=1; %Allows only a single value to be selected.
        end   
        SubtractMaps{i,1}=imresize3(SubtractMaps{i,1},Resample);
        SubtractMaps{i,2}=fileN;
        ParcelInfo.(PassStr).SubtractMaps{i,1}=fileN;
        if strcmpi(flipMaps,'yes')
            SubtractMaps{i,1}=flip(SubtractMaps{i,1},1);
        end      
        ParcelInfo.(PassStr).SubtractMaps{i,2}=SubtractMaps{i,1};
    end
else
    SubtractMaps=[];
    ParcelInfo.(PassStr).SubtractMaps=[];
end

numIntersectMaps=uiEnterName('0',['Enter # of intersection-maps',newline,'to load']);
numIntersectMaps=str2num(numIntersectMaps);
if numIntersectMaps~=0
    IntersectMaps=cell(numIntersectMaps,3);
    ParcelInfo.(PassStr).IntersectMaps=cell(numIntersectMaps,3);
    for i = 1:numIntersectMaps
        [fileN,dirN]=uigetfile([dirN,'*.*'],'Select intersection map (.nii .brik or .mat)');
        if contains(fileN,'.nii','IgnoreCase',true)
            intMap=load_nii([dirN,fileN]);
            IntersectMaps{i,1}=intMap.img;
        elseif contains(fileN,'.brik','IgnoreCase',true)
            [~,IntersectMaps{i,1},info] = BrikLoad ([dirN,fileN]);
        elseif contains(fileN,'.mat','IgnoreCase',true)
            SingleSelect=1; %Allows only a single value to be selected.
        end 
        IntersectMaps{i,1}=imresize3(IntersectMaps{i,1},Resample);
        IntersectMaps{i,2}=fileN;
        ParcelInfo.(PassStr).IntersectMaps{i,1}=fileN;
        if strcmpi(flipMaps,'yes')
            IntersectMaps{i,1}=flip(IntersectMaps{i,1},1);
        end        
        ParcelInfo.(PassStr).IntersectMaps{i,2}=IntersectMaps{i,1};
    end
else
    IntersectMaps=[];
    ParcelInfo.(PassStr).IntersectMaps=[];
end

if strcmpi(editParcel,'no')
    count = 1;
    UseMask=brainMap*0;
    UseLabels=cell(1);
end

ParcelNames=uiNameSelect(ParcelNames,'Select reference parcellation for ROI labels.',1);
if strcmpi(BrainSpace,'mni')
    load('Parcellations/anat/MNI152_1mmMNI.mat'); %mni_Anat;
    anatMap=imresize3(mni_Anat,Resample);
elseif strcmpi(BrainSpace,'tal')
    load('Parcellations/anat/TT_N27_1mmTAL.mat'); %tal_Anat;
    anatMap=imresize3(tal_Anat,Resample);
end

RefMap=load(['Parcellations/',ParcelNames]);

setClusterThresh=0;
while setClusterThresh==0  
    setBrainThresh=0; 
    brainMap=OrigBrainMap;
    [brainCoords,brainVals] = mat2coords(brainMap);       
    figure
    histogram(brainVals);
    hold on
    title('Distribution of brain values');
    
    while setBrainThresh==0

        setClusterThresh=0;    
        threshold=uiEnterName('1.96',['Select voxel-wise threshold',newline,'(use histogram for reference)']);
        threshold=str2num(threshold);
        tempbrainMap=single(brainMap>=threshold);
        figure
        imshow3D(tempbrainMap*4+anatMap);
        [KeepBrainThresh] = uiNameSelect({'Good','Try another threshold'},'Keep or change threshold:',1);
        if strcmpi(KeepBrainThresh,'Good')
            setBrainThresh=1; 
            brainMap=tempbrainMap;
            ParcelInfo.(PassStr).InitialMap{1,3}=threshold;
        end           
        close
    end
    close
    
    if ~isempty(SubtractMaps)
        setSubThresh=0;
        while setSubThresh==0
            tempSubtractMap=brainMap*0;
            SubThresholds=cell(numSubtractMaps,1);
            tempbrainMap=brainMap;
            for i = 1:numSubtractMaps
                [~,subVals] = mat2coords(SubtractMaps{i,1}); 
                figure
                histogram(subVals);
                hold on
                title(['Distribution of subtraction map values for:',newline,SubtractMaps{i,2}]);
                txtSubThresholds=uiEnterName('1.96',['Select voxel-wise threshold',newline,'(use histogram for reference)']);
                SubThresholds{i,1}=str2num(txtSubThresholds);    
                tempSubtractMap=tempSubtractMap+single(SubtractMaps{i,1}>=SubThresholds{i,1});
                close
            end
            tempbrainMap(tempSubtractMap~=0)=0;
            figure
            imshow3D(tempbrainMap*4+anatMap);  
            [KeepSubThresh] = uiNameSelect({'Good','Try another threshold'},'Keep or change threshold:',1);
            if strcmpi(KeepSubThresh,'Good')
                setSubThresh=1; 
                brainMap=tempbrainMap;
                ParcelInfo.(PassStr).SubtractMaps(:,3)=SubThresholds;
            end 
            close
        end
    end
    
    if ~isempty(IntersectMaps)
        setIntThresh=0;
        while setIntThresh==0
            tempIntersectMap=brainMap*0;
            IntThresholds=cell(numIntersectMaps,1);
            tempbrainMap=brainMap;
            for i = 1:numIntersectMaps
                [~,IntVals] = mat2coords(IntersectMaps{i,1}); 
                figure
                histogram(IntVals);
                hold on
                title(['Distribution of intersection map values for:',newline,IntersectMaps{i,2}]);
                txtIntThresholds=uiEnterName('1.96',['Select voxel-wise threshold',newline,'(use histogram for reference)']);
                IntThresholds{i,1}=str2num(txtIntThresholds);    
                tempIntersectMap=tempIntersectMap+single(IntersectMaps{i,1}>=IntThresholds{i,1});
                close
            end
            tempbrainMap(tempIntersectMap==0)=0;
            figure
            imshow3D(tempbrainMap*4+anatMap);  
            [KeepIntThresh] = uiNameSelect({'Good','Try another threshold'},'Keep or change threshold:',1);
            if strcmpi(KeepIntThresh,'Good')
                setIntThresh=1; 
                brainMap=tempbrainMap;
                ParcelInfo.(PassStr).IntersectMaps(:,3)=IntThresholds;                
            end 
            close
        end
    end    
    [brainCoords,brainVals] = mat2coords(brainMap); 
    clustersizethreshold=uiEnterName('10',['Select cluster size threshold',newline,'(use histogram for reference)']);
    clustersizethreshold=str2num(clustersizethreshold);    
    [clusterdef] = uiNameSelect({'face','edge','point'},'Select cluster-belonging criterion:',1);     
    [Clusters]=ClusterIdentifierSimple(brainVals,brainCoords,'val',0.9,clustersizethreshold,clusterdef);
    ClusterMap=brainMap*0;
    
    for i = 1:size(Clusters,1)
        
        cMap=coords2mat(Clusters{i,1},ClusterMap*0,ones(Clusters{i,3},1));
        ClusterMap(cMap==1)=i;
    end
    
    tempClusterMap=reshape(scale01(ClusterMap(:))*5,size(anatMap))+2;
    figure
    imshow3D(tempClusterMap+anatMap);
    hold on
    title('Current Cluster Map')    
    
    [KeepClusters] = uiNameSelect({'Good','Try another threshold'},'Perform analysis by subject or by run:',1);
    if strcmpi(KeepClusters,'Good')
        setClusterThresh=1;        
    end  
    ParcelInfo.(PassStr).ClusterInfo.ClusterThresh=clustersizethreshold;
    ParcelInfo.(PassStr).ClusterInfo.ClusterDef=clusterdef;
    close all
end



RefMap.UseMask=imresize3(RefMap.UseMask,size(brainMap),'nearest');
IgnoredMask=brainMap*0;
IgnoredLabels=cell(1);
IgnoredCount=1;
for i = 1:size(Clusters,1)
    
    [Label1,Label2] = NearestClusterLabel(Clusters{i,1},RefMap.UseMask,RefMap.UseLabels);
    tempMap=coords2mat(Clusters{i,1},brainMap*0,ones(size(Clusters{i,1},1),1));    
    figure
    imshow3D(tempMap*5+anatMap);
    hold on
    title(['Inspect current cluster',newline,Label1,newline,Label2]);    
    tempLabel=uiEnterName([Label1,'_or_',Label2],['Enter ROI label',newline,'Leave blank to remove cluster']);
    if strcmpi(tempLabel,'')
        close all
        IgnoredMask(tempMap==1)=IgnoredCount;
        IgnoredLabels{IgnoredCount,1}=[Label1,'_or_',Label2];
        IgnoredCount=IgnoredCount+1;        
        continue
    end
    if strcmpi(voxelOverwrite,'no')
        tempMap(UseMask~=0)=0;
    end
    UseMask(tempMap==1)=count;
    UseLabels{count,1}=tempLabel;
    count=count+1;
end
ParcelInfo.(PassStr).UseMask=UseMask;
ParcelInfo.(PassStr).UseLabels=UseLabels;
ParcelInfo.(PassStr).IgnoredMask=IgnoredMask;
ParcelInfo.(PassStr).IgnoredLabels=IgnoredLabels;

if strcmpi(editParcel,'yes')
    SaveName=LoadParcelName;
else
    SaveName=uiEnterName('',['Enter parcellation save name']);
end
save(['Parcellations/',SaveName],'UseMask','UseLabels','ParcelInfo');

%brainMap=imresize3(brainMap,Resample);



end
