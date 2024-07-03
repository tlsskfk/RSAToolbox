function [NumStr] = num2str4filename(Num,DecPlace)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    DecPlace=[];
end
if ~isempty(DecPlace)
    Num=round(Num,DecPlace);
end
NumStr=num2str(Num);
NumStr=strrep(NumStr,'.','p');
end

