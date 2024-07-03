function [outCell] = strrepCell(inCell,oldStr,newStr)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
warning('off','MATLAB:strrep:InvalidInputType');
if istable(inCell)
    inTable=inCell;
    inCell=table2cell(inTable);
    outCell=cellfun(@strrep,inCell,repmat({oldStr},size(inCell)),repmat({newStr},size(inCell)),'UniformOutput',0);
    inTable(:,:)=outCell;
    outCell=inTable;
else
    outCell=cellfun(@strrep,inCell,repmat({oldStr},size(inCell)),repmat({newStr},size(inCell)),'UniformOutput',0);
end

end

