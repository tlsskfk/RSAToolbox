function [OutNames] = uiEnterName(default_text,titletext)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin==0
    titletext='Enter name below';
    default_text='';
end
if nargin==1
    titletext='Enter name below';
end
UIFig=uifigure('Position',[500,200,300,150],'Visible','On');
Title=uilabel(UIFig,'Position',[0,150-45,300,50],'Text',titletext,'HorizontalAlignment','Center');
Text=uitextarea(UIFig,'Position',[25,150-80,250,25],'Value',default_text,'HorizontalAlignment','Center','Visible','On');
EndButton=uibutton(UIFig,'Position',[100,10,100,40],'Text','Enter','ButtonPushedFcn',@(EndButton,event) NextButtonPush(UIFig) );
waitfor(UIFig,'Visible','Off');
OutNames=Text.Value{1,1};
close all force;
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end