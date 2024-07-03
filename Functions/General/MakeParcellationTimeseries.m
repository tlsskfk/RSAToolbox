function [timeSeriesbyParcel,timeSeriesbyVoxel] = MakeParcellationTimeseries( Mat,Parcellation,varargin)
%Written by David Rothlein
%Input:

%Mat = 4D timeseries mat where the 4 diminsions are [X,Y,Z,time] and values
%are BOLD signal

%Parcellation = 3D matrix of integer values where each non-zero integer indicates
%which ROI/node a given voxel belongs to. 0s, NANs and INFs are excluded.
%if parcellation = [], Only timeSeriesbyVoxel is computed.

%varargin: 
%format- MakeParcellationTimeseries(... ,'name1',value1,'name2',value2,...)

    %%name: 'Voxel_ModeRemove', value: flag where any value other than [] calls process
    %%to remove voxels based on mode frequency. Logic is that voxels that do not have
    %%reliable signal will have many repeated values. Helpful when
    %%normalization applied such that voxels w/o signal may be non-zero.
    %%Default is [].

    %%name: 'Voxel_ModeThresh', value: scalar between 0 and 1 indicating
    %%within each voxel, the proportion of volumes that have identical values 
    %%over total number of volumes (i.e. mode frequency/num. volumes).
    %%Values greater than threshold will be removed.
    %%Default is 0.5 meaning if 50% of timpepoints in a given voxel have identical values, 
    %%that voxel will be removed.

    %%name: 'Voxel_CV_filter', value: scalar threshold whereby all voxels that have a
    %%larger CV value are removed. CV=coeffiecient of variation or STD /abs(mean)
    %%Default is to not remove any voxels based on CV.

    %%name: 'Voxel_CensorNAN', value: value in timeseries to be replaced by NAN
    %%NANs will subsequently be ignored or if specified, replaced using 'Voxel_Impute' 
    %%Default is 0. To skip this step, set value to []; 

    %%name: 'Voxel_Impute', value: flag where any value other than [] calls process
    %%to replace NANs with values using method indicated by 'Voxel_Impute_Type'
    %%Default is [].

    %%name: 'Voxel_Impute_Type', value: string indicating impute method options
    %%include-'previous','next','nearest','linear','spline','pchip'
    %%Default is 'linear'; see MATLAB function FILLMISSING for more info

    %%name: 'Voxel_Outlier', value: flag where any value other than [] calls process
    %%to identify and replace outliers with values at the outlier boundaries using 'clip' method.
    %%Default is []. see MATLAB function FILLOUTLIERS for more info

    %%name: 'Voxel_Outlier_Type', value: string indicating method to
    %%identify outliers
    %%options include-'median','mean','quartiles','grubbs','gesd'
    %%Default is 'mean' (mean indicates threshold by Z-score) 

    %%name: 'Voxel_Outlier_Thresh', value: scalar indicating outlier
    %%threshold value.
    %%Default is 3; for 'mean' = 3 = values greater or less than 3 standard
    %%deviations from the mean.

    %%name: 'VolumeIndex', value: 1 by L vector where L = number of volumes.
    %%Vector consists of 0 and 1s and 0s indicate volumes to be removed. 
    %%IMPORTANT: volumes are removed AFTER other steps to clean timeseries.
    %%Default is to not remove any volumes.

[ VolumeIndex ] = VariableSetter( 'VolumeIndex',[ones(1,size(Mat,4))],varargin);
[ TimeSeriesParams ] = VariableSetter( 'TimeSeriesParams',[],varargin); 
[ Voxel_ModeRemove ] = VariableSetter( 'Voxel_ModeRemove',[],varargin); 
[ Voxel_ModeThresh ] = VariableSetter( 'Voxel_ModeThresh',0.5,varargin);
[ Voxel_CV_filter ] = VariableSetter( 'Voxel_CV_filter',[],varargin);
[ Voxel_CensorNAN ] = VariableSetter( 'Voxel_CensorNAN',0,varargin);
[ Voxel_Impute ] = VariableSetter( 'Voxel_Impute',[],varargin);
[ Voxel_Impute_Type ] = VariableSetter( 'Voxel_Impute_Type','linear',varargin);
[ Voxel_Outlier ] = VariableSetter( 'Voxel_Outlier',[],varargin);
[ Voxel_Outlier_Type ] = VariableSetter( 'Voxel_Outlier_Type','mean',varargin);
[ Voxel_Outlier_Thresh ] = VariableSetter( 'Voxel_Outlier_Thresh',3,varargin);
Mat=double(Mat);
Parcellation=double(Parcellation);
if ~isempty(TimeSeriesParams)
    Voxel_ModeRemove=TimeSeriesParams.Voxel_ModeRemove;
    Voxel_ModeThresh=TimeSeriesParams.Voxel_ModeThresh;
    Voxel_CV_filter=TimeSeriesParams.Voxel_CV_filter;
    Voxel_CensorNAN=TimeSeriesParams.Voxel_CensorNAN;
    Voxel_Impute=TimeSeriesParams.Voxel_Impute;
    Voxel_Impute_Type=TimeSeriesParams.Voxel_Impute_Type;
    Voxel_Outlier=TimeSeriesParams.Voxel_Outlier;
    Voxel_Outlier_Type=TimeSeriesParams.Voxel_Outlier_Type;
    Voxel_Outlier_Thresh=TimeSeriesParams.Voxel_Outlier_Thresh;
end

Mat(isinf(Mat))=nan;
xL=size(Mat,1);
yL=size(Mat,2);
zL=size(Mat,3);
numVol=size(Mat,4);
skipParcel=0;
if isempty(Parcellation)
    Parcellation=ones(xL,yL,zL);
    skipParcel=1;
    timeSeriesbyParcel=[];
end
Parcellation(isnan(Parcellation))=0;
Parcellation(isinf(Parcellation))=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Step 1 Voxel Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mode filter = Remove voxels that do not vary  [Optional]
if ~isempty(Voxel_ModeRemove)
    [~,ModeFilter]=mode(abs(Mat),4);
    ModeFilter=ModeFilter/numVol;
    ModeFilter=single(ModeFilter<Voxel_ModeThresh);
else
    ModeFilter=ones(xL,yL,zL);
end

%CVfilter = Remove voxels that vary too much [Optional]
if ~isempty(Voxel_CV_filter)
    CVFilter=single((nanstd(Mat,0,4)./abs(nanmean(Mat,4)))<Voxel_CV_filter);
else
    CVFilter=ones(xL,yL,zL);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2: Reshape matricies from 4D to 2D. Match 2D timeseries mat %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove voxels w/ no signal [Not optional]
ZeroFilter=single(nanmean(abs(Mat),4)~=0);

%Select coordinates of usable voxels
FilterMat=single((ZeroFilter+CVFilter+ModeFilter)==3);
coords=mat2coords(FilterMat);

%Reshape TimeSeries mat and parcellation mat, ensuring the 1st dimension
%is aligned with coords
FilterMat=coords2mat(coords,FilterMat,[1:size(coords,1)]');
[coords,ind]=mat2coords(FilterMat);
FilterMat=reshape(FilterMat,[xL*yL*zL,1]);
timeSeriesMat=reshape(Mat,[xL*yL*zL,numVol]);
Parcellation=reshape(Parcellation,[xL*yL*zL,1]);
timeSeriesMat(FilterMat==0,:)=[];
Parcellation(FilterMat==0,:)=[];
FilterMat(FilterMat==0,:)=[];
if sum(single(ind==FilterMat))~=length(ind)
    disp('Coord mismatch. Inspect data')
end

timeSeriesMat(Parcellation==0,:)=[];
coords(Parcellation==0,:)=[];
Parcellation(Parcellation==0,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Step 3: Clean timeseries %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Identify censored time points (Default = 0) as NANs
if ~isempty(Voxel_CensorNAN)
    timeSeriesMat(timeSeriesMat==Voxel_CensorNAN)=nan;
end

%Impute NAN values
if ~isempty(Voxel_Impute)
    timeSeriesMat=fillmissing(timeSeriesMat,Voxel_Impute_Type,2);
else
    RemoveTimePts=sum(single(isnan(timeSeriesMat)))==size(timeSeriesMat,1);
    VolumeIndex(RemoveTimePts)=0;
end
%Replace Outliers
if ~isempty(Voxel_Outlier)
    timeSeriesMat=filloutliers(timeSeriesMat,'clip',Voxel_Outlier_Type,2,'ThresholdFactor',Voxel_Outlier_Thresh);
end

%Select a subset of volumes (defaults selects all volumes)

timeSeriesMat(:,VolumeIndex==0)=[];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Step 4: Create timeseries %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create timeseries by Voxel structure
timeSeriesbyVoxel.Coords=coords;
timeSeriesbyVoxel.Parcellation=Parcellation;
timeSeriesbyVoxel.Timeseries=timeSeriesMat;
ROICoords=[];
%Create timeseries by parcellation ROI/node structure
if skipParcel==0
    numROIs=max(Parcellation(:));
    timeSeriesParcelMat=zeros(numROIs,size(timeSeriesMat,2));
    count=1;
    for i = unique(Parcellation)'

        tempTimeSeries=timeSeriesMat(Parcellation==i,:);
        tempCoords=single(coords(Parcellation==i,:));
        ROICoords.CenterOfMass(count,:)=mean(tempCoords,1);
        ROICoords.STdev(count,:)=std(tempCoords,0,1);
        ROICoords.min(count,:)=min(tempCoords,[],1);
        ROICoords.max(count,:)=max(tempCoords,[],1);
        ROICoords.size(count,1)=size(tempCoords,1);

        timeSeriesParcelMat(count,:)=nanmean(tempTimeSeries,1);
        count=count+1;
    end
    timeSeriesParcelMat(isinf(timeSeriesParcelMat))=nan;
    if ~isempty(Voxel_Impute)
        timeSeriesParcelMat=fillmissing(timeSeriesParcelMat,Voxel_Impute_Type,2);
    end
    timeSeriesbyParcel.Timeseries=timeSeriesParcelMat;
    
    timeSeriesbyParcel.ROIinfo=ROICoords;
    timeSeriesbyParcel.Parcellation=unique(Parcellation);
end

