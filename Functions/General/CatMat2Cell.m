function [OutCell] = CatMat2Cell(CatMat,CatInd)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numCells=size(CatInd,1);
OutCell=cell(numCells,1);
for i = 1:numCells
    OutCell{i,1}=CatMat(CatInd(i,1):CatInd(i,2),:);
end
end

