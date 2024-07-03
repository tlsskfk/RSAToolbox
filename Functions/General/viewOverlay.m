function viewOverlay(ViewMap,RawVals)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    RawVals=0;
end

load('Parcellations/anat/MNI152_1mmMNI.mat'); %mni_Anat;
anatMap=imresize3(mni_Anat,size(ViewMap));
if RawVals==0
    tempClusterMap=reshape(scale01(ViewMap(:))*5,size(anatMap))+2;
end
figure
imshow3D(tempClusterMap+anatMap);

end

