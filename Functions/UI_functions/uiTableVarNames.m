function [NewVarOrder,NewVarNames] = uiTableVarNames(VarNames)
%Written by David Rothlein
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
VarNames=VarNames(:)';
UIFig=uifigure('Position',[500,200,600,425],'Visible','On');
Title=uilabel(UIFig,'Position',[20,400,300,30],'Text','Edit table variable order and names','HorizontalAlignment','Left');
ItemCounter=uislider(UIFig,'Visible','Off','Value',1);
NewVarNameText=uitextarea(UIFig,'Position',[150,60,300,25],'Value','Edit new variable name here','HorizontalAlignment','Center','Visible','On');    
numItems=length(VarNames);
EmptyItems=cell(1,numItems);
for i = 1:numItems
    EmptyItems{1,i}='';
end
VarNamesListBox=uilistbox(UIFig,'Position',[50,100,200,250],'Items',VarNames,'Multiselect','Off','Visible','On');    
VarLabel=uilabel(UIFig,'Position',[50,355,200,30],'Text','Select variable name to add and edit','HorizontalAlignment','Center');

VarOrderListBox=uilistbox(UIFig,'Position',[350,250,200,100],'Items',EmptyItems,'Visible','On');
NewVarNameListBox=uilistbox(UIFig,'Position',[350,100,200,100],'Items',EmptyItems,'Visible','On');  
VarOrderLabel=uilabel(UIFig,'Position',[350,355,200,30],'Text','New Order','HorizontalAlignment','Center');
VarNameLabel=uilabel(UIFig,'Position',[350,205,200,30],'Text','New Order and Name','HorizontalAlignment','Center');

EndButton=uibutton(UIFig,'Position',[450,10,100,40],'Text','End','Visible','Off','ButtonPushedFcn',@(EndButton,event) EndButtonPush(UIFig) );
UndoButton=uibutton(UIFig,'Position',[50,10,100,40],'Text','Undo','Visible','Off','ButtonPushedFcn',@(UndoButton,event) UndoButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameListBox,NewVarNameText,ItemCounter,EndButton) );
AddNewNameButton=uibutton(UIFig,'Position',[270,150,60,30],'Visible','Off','Text','Add Name->','ButtonPushedFcn',@(AddNewNameButton,event) NewNameButtonPush(NewVarNameText,NewVarNameListBox,AddNewNameButton,ItemCounter,VarNamesListBox,EndButton,UndoButton) );
AddNewOrderButton=uibutton(UIFig,'Position',[270,300,60,30],'Text','Add->','ButtonPushedFcn',@(AddNewOrderButton,event) VarOrderButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameText,AddNewNameButton,ItemCounter,UndoButton) );
UseDefaultButton=uibutton(UIFig,'Position',[260,10,80,30],'Text','Use Default','ButtonPushedFcn',@(UseDefaultButton,event) UseDefaultButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameListBox,ItemCounter,UndoButton,EndButton) );

waitfor(UIFig,'Visible','Off');
NewVarOrder=VarOrderListBox.Items;
NewVarNames=NewVarNameListBox.Items;
close all force;


end

function VarOrderButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameText,AddNewNameButton,ItemCounter,UndoButton)
    VarOrderListBox.Items{ItemCounter.Value}=VarNamesListBox.Value;
    ItemCounter.Value=ItemCounter.Value+1;
    NewVarNameText.Value=VarNamesListBox.Value;
    NewVarNameText.Value=strrep(NewVarNameText.Value,'_','');
    NewVarNameText.Value=strrep(NewVarNameText.Value,'.mat','');
    VarCounter=find(ismember(VarNamesListBox.Items,VarNamesListBox.Value));
    VarNamesListBox.Items(VarCounter)=[];
    AddNewNameButton.Visible='On';
    UndoButton.Visible='Off';
end

function NewNameButtonPush(NewVarNameText,NewVarNameListBox,AddNewNameButton,ItemCounter,VarNamesListBox,EndButton,UndoButton)
    NewVarNameListBox.Items{1,ItemCounter.Value-1}=NewVarNameText.Value{1,1};
    NewVarNameText.Value='';
    AddNewNameButton.Visible='Off';
    UndoButton.Visible='On';
    if isempty(VarNamesListBox.Items)
        EndButton.Visible='On';
    end
end

function UndoButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameListBox,NewVarNameText,ItemCounter,EndButton,AddNewOrderButton)
    if ItemCounter.Value > 1
        ItemCounter.Value=ItemCounter.Value-1;
        VarNamesListBox.Items=[VarNamesListBox.Items,VarOrderListBox.Items(ItemCounter.Value)];
        VarOrderListBox.Items{ItemCounter.Value}='';
        NewVarNameListBox.Items{ItemCounter.Value}='';        
        NewVarNameText.Value='';
        EndButton.Visible='Off';
        AddNewOrderButton.Visible='On';
    end
end

function UseDefaultButtonPush(VarNamesListBox,VarOrderListBox,NewVarNameListBox,ItemCounter,UndoButton,EndButton)
    VarOrderListBox.Items=VarNamesListBox.Items;
    NewVarNameListBox.Items=VarNamesListBox.Items;
    NewVarNameListBox.Items=strrepCell(NewVarNameListBox.Items,'.mat','');
    EndButton.Visible='On';
    UndoButton.Visible='On';
end

function EndButtonPush(UIFig)
    UIFig.Visible='Off';
end