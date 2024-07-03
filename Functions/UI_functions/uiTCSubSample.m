function [SubsampleParams,SaveString] = uiTCSubSample(SplitByNames)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
UIFig=uifigure('Position',[500,50,300,500],'Visible','On');
Title=uilabel(UIFig,'Position',[0,500-45,300,50],'Text','Set Subsample Parameters','HorizontalAlignment','Center');
EndButton=uibutton(UIFig,'Position',[100,10,100,40],'Text','Done','ButtonPushedFcn',@(EndButton,event) NextButtonPush(UIFig) );

NumSubsTitle=uilabel(UIFig,'Position',[100,500-70,130,30],'Text',['Enter # (odd integer)',newline,'of subsamples'],'HorizontalAlignment','left','Visible','On');
NumSubsText=uitextarea(UIFig,'Position',[20,500-70,50,30],'Value','1','HorizontalAlignment','Center','Visible','On');

SubSampleSizeText=uitextarea(UIFig,'Position',[20,500-120,50,30],'Value','0.5','HorizontalAlignment','Center','Visible','On');
SubSampleSizeTitle=uilabel(UIFig,'Position',[100,500-120,130,30],'Text',['Enter subsample size',newline,'from 0 to 1'],'HorizontalAlignment','left');

SaveAppendCheck=uicheckbox(UIFig,'Position',[20,500-170,150,30],'Text','');
SaveAppendTitle=uilabel(UIFig,'Position',[100,500-170,130,30],'Text',['Append to saved',newline,'RSMs?'],'HorizontalAlignment','left');

NumRepsText=uitextarea(UIFig,'Position',[20,500-220,50,30],'Value','1','HorizontalAlignment','Center','Visible','On');
NumRepsTitle=uilabel(UIFig,'Position',[100,500-220,200,30],'Text',['Enter # (integer)',newline,'of repetitions'],'HorizontalAlignment','left');

SplitTCTitle=uilabel(UIFig,'Position',[50,500-260,200,20],'Text','Select split variable:','HorizontalAlignment','center');
SplitTCListBox=uilistbox(UIFig,'Position',[50,500-430,200,170],'Items',SplitByNames,'Multiselect','off','Value',{});

waitfor(UIFig,'Visible','Off');
%%%% Set param values here %%%%

SubsampleParams.NumSubs=str2num(NumSubsText.Value{1,1});
SubsampleParams.SaveAppend=single(SaveAppendCheck.Value);
SubsampleParams.SplitVar=SplitTCListBox.Value;
SubsampleParams.NumReps=str2num(NumRepsText.Value{1,1});
SubsampleParams.SubSampleSize=str2num(SubSampleSizeText.Value{1,1});
if mod(SubsampleParams.NumSubs,2)~=1
    SubsampleParams.NumSubs=SubsampleParams.NumSubs+1;
end
SaveString=['TC-',SubsampleParams.SplitVar,'_numSubs-',NumSubsText.Value{1,1},'_subSize-',num2str4filename(SubsampleParams.SubSampleSize,2)];
end


function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end