function [ROIs,Coords] = uiROICoordAssigment(VarNames,TitleText)
%Written by David Rothlein
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    TitleText = 'Edit table variable order and names';
end
VarNames=VarNames(:)';
UIFig=uifigure('Position',[500,200,800,425],'Visible','On');
Title=uilabel(UIFig,'Position',[20,400,300,30],'Text',TitleText,'HorizontalAlignment','Left');
ItemCounter=uislider(UIFig,'Visible','Off','Value',0);
XLabel=uilabel(UIFig,'Position',[50,75,50,25],'Text','X','HorizontalAlignment','Center');
YLabel=uilabel(UIFig,'Position',[125,75,50,25],'Text','Y','HorizontalAlignment','Center');
ZLabel=uilabel(UIFig,'Position',[200,75,50,25],'Text','Z','HorizontalAlignment','Center');
XCoordNameText=uitextarea(UIFig,'Position',[50,50,50,25],'Value','','HorizontalAlignment','Center','Visible','On');    
YCoordNameText=uitextarea(UIFig,'Position',[125,50,50,25],'Value','','HorizontalAlignment','Center','Visible','On');    
ZCoordNameText=uitextarea(UIFig,'Position',[200,50,50,25],'Value','','HorizontalAlignment','Center','Visible','On');    
numItems=length(VarNames);
EmptyItems=cell(1,numItems);
count=1;
Data={'','','',''};
colNames={'ROI','X','Y','Z'};
for i = 1:numItems
    EmptyItems{1,i}='';
end
VarNamesListBox=uilistbox(UIFig,'Position',[50,100,200,250],'Items',VarNames,'Multiselect','Off','Visible','On');    
VarLabel=uilabel(UIFig,'Position',[50,355,200,30],'Text','Select variable name to add and edit','HorizontalAlignment','Center');

VarOrderListBox=uitable(UIFig,'Position',[350,100,400,250],'Data',Data,'ColumnName',colNames,'Visible','On');
VarOrderLabel=uilabel(UIFig,'Position',[350,355,200,30],'Text','ROIs and coords','HorizontalAlignment','Center');

EndButton=uibutton(UIFig,'Position',[600,10,100,40],'Text','End','Visible','Off','ButtonPushedFcn',@(EndButton,event) EndButtonPush(UIFig) );
UndoButton=uibutton(UIFig,'Position',[400,10,100,40],'Text','Undo','Visible','Off','ButtonPushedFcn',@(UndoButton,event) UndoButtonPush(VarNamesListBox,VarOrderListBox,ItemCounter,EndButton) );
AddNewNameButton=uibutton(UIFig,'Position',[270,50,100,40],'Visible','Off','Text','Add Coords->','ButtonPushedFcn',@(AddNewNameButton,event) NewNameButtonPush(VarOrderListBox,AddNewNameButton,ItemCounter,VarNamesListBox,EndButton,UndoButton,XCoordNameText,YCoordNameText,ZCoordNameText) );
AddNewOrderButton=uibutton(UIFig,'Position',[270,300,60,30],'Text','Add->','ButtonPushedFcn',@(AddNewOrderButton,event) VarOrderButtonPush(VarNamesListBox,VarOrderListBox,AddNewNameButton,ItemCounter,UndoButton) );

waitfor(UIFig,'Visible','Off');
AllData=VarOrderListBox.Data;
Coords=cell2mat(cellfun(@str2num, AllData(:,[2:4]),'UniformOutput',0));
ROIs=AllData(:,1);
close all force;


end

function VarOrderButtonPush(VarNamesListBox,VarOrderListBox,AddNewNameButton,ItemCounter,UndoButton)
    ItemCounter.Value=ItemCounter.Value+1;
    VarOrderListBox.Data{ItemCounter.Value,1}=VarNamesListBox.Value;
    VarCounter=find(ismember(VarNamesListBox.Items,VarNamesListBox.Value));
    VarNamesListBox.Items(VarCounter)=[];
    AddNewNameButton.Visible='On';
    UndoButton.Visible='Off';
end

function NewNameButtonPush(VarOrderListBox,AddNewNameButton,ItemCounter,VarNamesListBox,EndButton,UndoButton,XCoordNameText,YCoordNameText,ZCoordNameText)
    VarOrderListBox.Data{ItemCounter.Value,2}=XCoordNameText.Value{1,1};
    VarOrderListBox.Data{ItemCounter.Value,3}=YCoordNameText.Value{1,1};
    VarOrderListBox.Data{ItemCounter.Value,4}=ZCoordNameText.Value{1,1};
    XCoordNameText.Value{1,1}='';
    YCoordNameText.Value{1,1}='';
    ZCoordNameText.Value{1,1}='';
    AddNewNameButton.Visible='Off';
    UndoButton.Visible='On';
    if isempty(VarNamesListBox.Items)
        EndButton.Visible='On';
    end
end

function UndoButtonPush(VarNamesListBox,VarOrderListBox,ItemCounter,EndButton)
    if ItemCounter.Value > 0
        ItemCounter.Value=ItemCounter.Value;
        VarNamesListBox.Items=[VarNamesListBox.Items,VarOrderListBox.Data{ItemCounter.Value,1}];
        VarOrderListBox.Data(ItemCounter.Value,:)=[];
        ItemCounter.Value=ItemCounter.Value-1;
        EndButton.Visible='Off';
    end
end

function EndButtonPush(UIFig)
    UIFig.Visible='Off';
end