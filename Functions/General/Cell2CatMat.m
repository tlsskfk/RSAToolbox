function [OutMat,CatInd] = Cell2CatMat(InCell)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
numCells=size(InCell,1);
CatInd=zeros(numCells,2);
StartInd=1;
OutMat=[];
for i =1:numCells
    TempMat=InCell{i,1};
    cellLength=size(TempMat,1);
    CatInd(i,:)=[StartInd,StartInd+cellLength-1];
    StartInd=StartInd+cellLength;
    OutMat=[OutMat;TempMat];
end    
end

