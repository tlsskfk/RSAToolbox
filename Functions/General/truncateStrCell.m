function [OutCell] = truncateStrCell(InCell,truncLength)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin ==1
    truncLength=namelengthmax;
end
for i = 1:length(InCell)
    if length(InCell{i}) > truncLength
        removeN=length(InCell{i})-truncLength;
        t=randperm(length(InCell{i}))';
        t(removeN+1:end,:)=[];
        tempStr=InCell{i};
        tempStr(:,t)=[];
        InCell{i}=tempStr;
    end
end
OutCell=InCell;    
end

