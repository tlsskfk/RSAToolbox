function [SearchCoordVal,minInd,minVal] = nearestCoordInd(inputCoord,searchCoordSet,SearchCoordVals)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
inputCoord=double(inputCoord);
searchCoordSet=double(searchCoordSet);
inputCoord=repmat(inputCoord,[size(searchCoordSet,1),1]);
[minVal,minInd]=min(sum(abs(inputCoord - searchCoordSet),2));
SearchCoordVal=SearchCoordVals(minInd,1);
end

