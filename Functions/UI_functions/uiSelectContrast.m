function [cMat,cNames] = uiSelectContrast(CondNames)
%Written by David Rothlein
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Next=1;
numConds=length(CondNames);
cMat=[];
cNames=cell(1);
CondNames=CondNames(:)';
count=1;
while Next==1 
    UIFig=uifigure('Position',[500,200,600,425],'Visible','On');
    Title=uilabel(UIFig,'Position',[20,400,300,30],'Text','Set GLM contrasts','HorizontalAlignment','Left');
    
    ContrastName=uitextarea(UIFig,'Position',[150,60,300,25],'Value','Enter Contrast Name Here','HorizontalAlignment','Center','Visible','On');    
 
    CondNamesListBox=uilistbox(UIFig,'Position',[50,100,200,250],'Items',CondNames,'Multiselect','On','Visible','On');    
    CondLabel=uilabel(UIFig,'Position',[50,355,200,30],'Text','Select conditions to add','HorizontalAlignment','Center');

    PosListBox=uilistbox(UIFig,'Position',[350,250,200,100],'Items',{''},'Visible','On');
    NegListBox=uilistbox(UIFig,'Position',[350,100,200,100],'Items',{''},'Visible','On');  
    PosLabel=uilabel(UIFig,'Position',[350,355,200,30],'Text','Positive','HorizontalAlignment','Center');
    NegLabel=uilabel(UIFig,'Position',[350,205,200,30],'Text','Negative','HorizontalAlignment','Center');
    
    
    PosButton=uibutton(UIFig,'Position',[270,300,60,30],'Text','Add->','ButtonPushedFcn',@(PosButton,event) AddButtonPush(CondNamesListBox,PosListBox) );
    NegButton=uibutton(UIFig,'Position',[270,150,60,30],'Text','Add->','ButtonPushedFcn',@(NegButton,event) AddButtonPush(CondNamesListBox,NegListBox) );
    ClearButton=uibutton(UIFig,'Position',[50,10,100,40],'Text','Clear','ButtonPushedFcn',@(ClearButton,event) ClearButtonPush(PosListBox,NegListBox) );
    NextButton=uibutton(UIFig,'Position',[250,10,100,40],'Text','Next','ButtonPushedFcn',@(NextButton,event) NextButtonPush(UIFig) );
    EndButton=uibutton(UIFig,'Position',[450,10,100,40],'Text','End','ButtonPushedFcn',@(EndButton,event) EndButtonPush(UIFig,CondNamesListBox) );
    
    waitfor(UIFig,'Visible','Off');
    if strcmp(CondNamesListBox.Value{1,1},'end2212')
        Next=0;
    end
    posMat=single(ismember(CondNames,PosListBox.Items));
    negMat=single(ismember(CondNames,NegListBox.Items));
    cMat = [cMat;posMat*sum(negMat(:))-negMat*sum(posMat)];
    cNames{count,1}=ContrastName.Value{1,1};
    count=count+1;
    close all force;
end

end

function AddButtonPush(CondNamesListBox,AddListBox)
    AddListBox.Items=CondNamesListBox.Value;
end

function ClearButtonPush(PosListBox,NegListBox)
    PosListBox.Items={''};
    NegListBox.Items={''};
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end

function EndButtonPush(UIFig,CondNamesListBox)
    CondNamesListBox.Items={'end2212'};
    CondNamesListBox.Value={'end2212'};
    UIFig.Visible='Off';
end