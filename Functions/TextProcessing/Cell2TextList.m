function [TextList] = Cell2TextList(CellList)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
TextList=[];
CellList=CellList(:);
if length(CellList) == 1
    TextList=CellList{1,1};
elseif length(CellList) == 2
    TextList=[CellList{1,1},' and ',CellList{2,1}];
else
    for i = 1:length(CellList)
        if i == length(CellList)
            TextList=[TextList,' and ',CellList{i,1}];
        else
            TextList=[TextList,CellList{i,1},', '];
        end
    end
end

