function [OutTable] = WordStudy_pRMSConvert(InRSMStruct,flipsign)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    flipsign=0;
end
fNames=fieldnames(InRSMStruct);
for i = 1:length(fNames)
    if flipsign == 1
        InRSMStruct.(fNames{i})=mat2uppertriuvectormat(triRSMs2SymRSMs(InRSMStruct.(fNames{i})))*-1;
    else
        InRSMStruct.(fNames{i})=mat2uppertriuvectormat(triRSMs2SymRSMs(InRSMStruct.(fNames{i})));
    end
end 
OutTable=struct2table(InRSMStruct);
end