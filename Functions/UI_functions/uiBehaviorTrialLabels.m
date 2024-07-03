function [LabelAssignments,Labels,ItemGroups] = uiBehaviorTrialLabels(ItemNames,UseTitle)
%Written by David Rothlein
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    UseTitle='Assign a condition label to a set of items';
end

ItemNames=unique(ItemNames);
ItemNames=ItemNames(:)';
LabelAssignments=cell(1);
Labels={''};
ItemGroups=cell(1);
Next=1;
count=1;
while Next==1 
    UIFig=uifigure('Position',[500,200,1100,425],'Visible','On');
    FigTitle=uilabel(UIFig,'Position',[20,400,300,30],'Text',UseTitle,'HorizontalAlignment','Left');
    
    ConditionName=uitextarea(UIFig,'Position',[400,60,300,25],'Value','Enter Condition Name Here','HorizontalAlignment','Center','Visible','On');    
 
    ItemNamesListBox=uilistbox(UIFig,'Position',[50,100,300,250],'Items',ItemNames,'Multiselect','On','Visible','On');    
    ItemLabel=uilabel(UIFig,'Position',[50,355,300,30],'Text','Select item(s) to label','HorizontalAlignment','Center');

    SavedNamesListBox=uilistbox(UIFig,'Position',[750,100,300,250],'Items',Labels,'Multiselect','On','Visible','On');    
    SavedLabel=uilabel(UIFig,'Position',[750,355,300,30],'Text','Saved Labels','HorizontalAlignment','Center');
    
    ItemGroupListBox=uilistbox(UIFig,'Position',[450,100,200,250],'Items',{''},'Visible','On');
    ItemGroupLabel=uilabel(UIFig,'Position',[450,355,200,30],'Text','Items assigned to condition label','HorizontalAlignment','Center');
       
    AddButton=uibutton(UIFig,'Position',[370,200,60,30],'Text','Add->','ButtonPushedFcn',@(AddButton,event) AddButtonPush(ItemNamesListBox,ItemGroupListBox,ConditionName) );
    ClearButton=uibutton(UIFig,'Position',[150,10,100,40],'Text','Clear','ButtonPushedFcn',@(ClearButton,event) ClearButtonPush(ItemGroupListBox) );
    NextButton=uibutton(UIFig,'Position',[500,10,100,40],'Text','Next','ButtonPushedFcn',@(NextButton,event) NextButtonPush(UIFig) );
    EndButton=uibutton(UIFig,'Position',[850,10,100,40],'Text','End','ButtonPushedFcn',@(EndButton,event) EndButtonPush(UIFig,ItemNamesListBox) );
    DefaultButton=uibutton(UIFig,'Position',[650,10,100,40],'Text','Default','ButtonPushedFcn',@(DefaultButton,event) DefaultButtonPush(UIFig,FigTitle) );
    waitfor(UIFig,'Visible','Off');
    if strcmp(ItemNamesListBox.Value{1,1},'end2212')
        Next=0;
    end
    if strcmpi(FigTitle.Text,UseTitle)
        LabelAssignments{1,count}{1,1}=ConditionName.Value{1,1};
        LabelAssignments{1,count}{1,2}=ItemGroupListBox.Items(:)';
        Labels{1,count}=ConditionName.Value{1,1};
        ItemGroups{1,count}=ItemGroupListBox.Items(:)';
        
        count=count+1;
    else
        for i = 1:length(ItemNamesListBox.Items)
            LabelAssignments{1,i}{1,1}=ItemNamesListBox.Items{1,i};
            LabelAssignments{1,i}{1,2}=ItemNamesListBox.Items(1,i);
            Labels{1,i}=ItemNamesListBox.Items{1,i};
            ItemGroups{1,i}=ItemNamesListBox.Items(1,i);    
        end
        Next=0;
    end
    close all force;
end

end

function AddButtonPush(CondNamesListBox,AddListBox,ConditionName)
    AddListBox.Items=CondNamesListBox.Value;
    tempName=AddListBox.Items{1,1};
    tempDotInd=strfind(tempName,'.');
    if isempty(tempDotInd)
        ConditionName.Value=tempName;
    else
        ConditionName.Value=tempName(1,1:tempDotInd(end)-1);
    end
end

function ClearButtonPush(PosListBox)
    PosListBox.Items={''};
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end

function EndButtonPush(UIFig,CondNamesListBox)
    CondNamesListBox.Items={'end2212'};
    CondNamesListBox.Value={'end2212'};
    UIFig.Visible='Off';
end

function DefaultButtonPush(UIFig,FigTitle)
    FigTitle.Text='end2212';
    UIFig.Visible='Off';
end