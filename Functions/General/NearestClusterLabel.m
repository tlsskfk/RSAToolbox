function [Label1,Label2] = NearestClusterLabel(ROICoords,ReferenceMap,ReferenceLabels)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    ROICoords=single(ROICoords);
    AveCoords=mean(ROICoords,1);
    ReferenceLabels=ReferenceLabels(:);
    [RefCoords,RefVals]=mat2coords(ReferenceMap);
    [SearchCoordVal,~,~] = nearestCoordInd(AveCoords,RefCoords,RefVals);
    Label1=ReferenceLabels{SearchCoordVal,1};
    tempRefCoords=RefCoords(RefVals~=SearchCoordVal,:);
    tempRefVals=RefVals(RefVals~=SearchCoordVal,:);
    [SearchCoordVal,~,~] = nearestCoordInd(AveCoords,tempRefCoords,tempRefVals);
    Label2=ReferenceLabels{SearchCoordVal,1};

end

