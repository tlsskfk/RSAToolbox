function [f] = plotByRun(subs,plotData)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
subNames=unique(subs);
figure
for i =1:length(subNames)
    pData=plotData(ismember(subs,subNames{i,1}),1);
    f=plot([1:length(pData)]',pData);
    hold on
end   
end

