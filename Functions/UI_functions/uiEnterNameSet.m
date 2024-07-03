function [NameSet] = uiEnterNameSet(titleText)
%Written by David Rothlein
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin==0
    titleText='Enter names or labels';
end
NameSet=[];
UIFig=uifigure('Position',[500,200,600,425],'Visible','On');
Title=uilabel(UIFig,'Position',[20,400,300,30],'Text',titleText,'HorizontalAlignment','Left');
ItemCounter=uislider(UIFig,'Visible','Off','Value',1);
NewVarNameText=uitextarea(UIFig,'Position',[25,250,300,25],'Value','Enter name here','HorizontalAlignment','Center','Visible','On');    

VarNamesListBox=uilistbox(UIFig,'Position',[350,100,200,250],'Items',{''},'Visible','On');
VarNamesLabel=uilabel(UIFig,'Position',[350,355,200,30],'Text','Entered names','HorizontalAlignment','Center');

EndButton=uibutton(UIFig,'Position',[450,10,100,40],'Text','End','Visible','On','ButtonPushedFcn',@(EndButton,event) EndButtonPush(UIFig) );
UndoButton=uibutton(UIFig,'Position',[50,10,100,40],'Text','Undo','Visible','Off','ButtonPushedFcn',@(UndoButton,event) UndoButtonPush(VarNamesListBox,NewVarNameText,ItemCounter) );
AddNameButton=uibutton(UIFig,'Position',[135,175,80,50],'Text','Add->','ButtonPushedFcn',@(AddNameButton,event) AddNameButtonPush(VarNamesListBox,NewVarNameText,ItemCounter,UndoButton) );

waitfor(UIFig,'Visible','Off');
NameSet=VarNamesListBox.Items;
close all force;


end

function AddNameButtonPush(VarNamesListBox,NewVarNameText,ItemCounter,UndoButton)
    VarNamesListBox.Items{ItemCounter.Value}=NewVarNameText.Value{1,1};
    NewVarNameText.Value='';
    ItemCounter.Value=ItemCounter.Value+1;
    UndoButton.Visible='On';
end

function UndoButtonPush(VarNamesListBox,NewVarNameText,ItemCounter)
    if ItemCounter.Value > 1
        ItemCounter.Value=ItemCounter.Value-1;
        VarNamesListBox.Items(:,ItemCounter.Value)=[]; 
        NewVarNameText.Value='';
    end
end

function EndButtonPush(UIFig)
    UIFig.Visible='Off';
end