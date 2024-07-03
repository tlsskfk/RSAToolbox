function [SplitParams,SaveString] = uiTCSplit(SplitByNames)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
UIFig=uifigure('Position',[500,50,300,500],'Visible','On');
Title=uilabel(UIFig,'Position',[0,500-45,300,50],'Text','Set Timeseries Parameters','HorizontalAlignment','Center');
EndButton=uibutton(UIFig,'Position',[100,10,100,40],'Text','Done','ButtonPushedFcn',@(EndButton,event) NextButtonPush(UIFig) );

SaveRSMCheck=uicheckbox(UIFig,'Position',[20,500-70,150,30],'Text','');
SaveRSMCheckTitle=uilabel(UIFig,'Position',[100,500-70,130,30],'Text','Save Split RSMs?','HorizontalAlignment','left');

SaveSplitInfoCheck=uicheckbox(UIFig,'Position',[20,500-130,150,30],'Text','');
SaveSplitInfoCheckTitle=uilabel(UIFig,'Position',[100,500-130,130,30],'Text','Save Split Info?','HorizontalAlignment','left');

ComputeSplitRFCheck=uicheckbox(UIFig,'Position',[20,500-100,150,30],'Text','');
ComputeSplitRFTitle=uilabel(UIFig,'Position',[100,500-100,130,30],'Text','Compute Split RF','HorizontalAlignment','left');

NumRepsTitle=uilabel(UIFig,'Position',[100,500-430,200,30],'Text',['Enter number (integer) of',newline,'split iterations'],'HorizontalAlignment','left','Visible','On');
NumRepsText=uitextarea(UIFig,'Position',[20,500-430,50,30],'Value','1','HorizontalAlignment','Center','Visible','On');

SplitTCTitle=uilabel(UIFig,'Position',[50,500-180,200,20],'Text','Select split variable:','HorizontalAlignment','center');
SplitTCListBox=uilistbox(UIFig,'Position',[50,500-350,200,170],'Items',SplitByNames,'Multiselect','off','Value',{},'ValueChangedFcn',@(SplitTCListBox,event) SplitTCSelectFcn(SplitTCListBox.Value,NumRepsTitle,NumRepsText) );

waitfor(UIFig,'Visible','Off');
%%%% Set param values here %%%%

SplitParams.SaveRSM=single(SaveRSMCheck.Value);
SplitParams.SaveSplitInfo=single(SaveSplitInfoCheck.Value);
SplitParams.ComputeSplitRF=single(ComputeSplitRFCheck.Value);
SplitParams.SplitVar=SplitTCListBox.Value;
SplitParams.NumReps=str2num(NumRepsText.Value{1,1});

SaveString=['numReps-',NumRepsText.Value{1,1}];
end


function SplitTCSelectFcn(SplitName,NumRepsTitle,NumRepsText)
    if strcmpi(SplitName,'random')
        NumRepsTitle.Visible='On';
        NumRepsText.Visible='On';
    end    
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end