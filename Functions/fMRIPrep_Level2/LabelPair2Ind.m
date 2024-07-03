function [ ind,reverseInd ] = LabelPair2Ind( Labels,LabelPairs,delim )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    delim='_';
end
dL=length(delim);
ind=cell(length(LabelPairs),1);
reverseInd=ind;
for i = 1:length(LabelPairs)
   lPair=LabelPairs{i,1};
   x=strfind(lPair,delim);
   a=lPair(1:x-1);
   b=lPair(x+dL:end);
   ind1=find(ismember(Labels,a));
   ind2=find(ismember(Labels,b));
   if isempty(ind2)
       ind2=find(contains(Labels,b));
   end
   ind{i,1}=[ind1,ind2];
   reverseInd{i,1}=[ind2,ind1];
end

