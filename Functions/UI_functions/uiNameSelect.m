function [OutNames] = uiNameSelect(inNames,titletext,singleSelect,noForceClose)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    titletext='Select/highlight names below and click next';
    singleSelect=0;
    noForceClose=0;
end
if nargin==2
    singleSelect=0;
    noForceClose=0;
end
if nargin==3
    noForceClose=0;
end
if length(inNames)>1
    UIFig=uifigure('Position',[500,50,300,500],'Visible','On');
    Title=uilabel(UIFig,'Position',[0,500-45,300,50],'Text',titletext,'HorizontalAlignment','Center');
    if singleSelect==0
        ListBox=uilistbox(UIFig,'Position',[20,60,260,400],'Items',inNames,'Multiselect','on','Value',{});
    else
        ListBox=uilistbox(UIFig,'Position',[20,60,260,400],'Items',inNames,'Multiselect','off','Value',{});
    end
    EndButton=uibutton(UIFig,'Position',[100,10,100,40],'Text','Next','ButtonPushedFcn',@(EndButton,event) NextButtonPush(UIFig) );
    waitfor(UIFig,'Visible','Off');
    OutNames=ListBox.Value;
    if iscell(OutNames)
        OutNames=OutNames';
    end  
    if noForceClose==0
        close all force;
    end
else
    OutNames=inNames;
end
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end