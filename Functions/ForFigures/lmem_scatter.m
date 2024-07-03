function [fs] = lmem_scatter(groupingVar,y,x,y_label,x_label,ScatterTitle,statsText)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    y_label='';
    x_label='';
    ScatterTitle='';
    statsText='';
end
if nargin==4
    x_label='';    
    statsText='';
    ScatterTitle='';    
end
if nargin==5  
    statsText='';
    ScatterTitle='';    
end
if nargin==6  
    statsText='';   
end

if ~isempty(groupingVar)
    figure;
    Groups=unique(groupingVar);
    numGroups=length(Groups);
    fs=cell(2,numGroups);
    for i = 1:numGroups
        selectInd=ismember(groupingVar,Groups{i,1});
        t_y=y(selectInd,1);
        t_x=x(selectInd,1);
        b=regress(t_x,[t_y,ones(length(t_y),1)]);
        slope=b(1,1);
        Int=b(2,1);
        xRange=[min(t_y),max(t_y)];
        yRange=[Int+xRange(1,1)*slope,Int+xRange(1,2)*slope];
        fs{1,i}=scatter(t_y,t_x,'filled');
        hold on
        fs{2,i}=plot(xRange,yRange,'LineWidth',2);
        fs{2,i}.Color=fs{1,i}.CData;
        hold on
    end
    b=regress(x,[y,ones(length(y),1)]);
    slope=b(1,1);
    Int=b(2,1);
    xRange=[min(y),max(y)];
    yRange=[Int+xRange(1,1)*slope,Int+xRange(1,2)*slope];
    plot(xRange,yRange,'k','LineWidth',3);
    title(ScatterTitle,'Interpreter','none');
    ylims=ylim;
    xlims=xlim;
    yUnit=(ylims(1,2)-ylims(1,1))/100;
    xUnit=(xlims(1,2)-xlims(1,1))/100;
    if slope > 0
        text(xUnit*90,ylims(1,1)+yUnit*15,statsText,'HorizontalAlignment','Center','Interpreter','none','FontSize',8);
    else
        text(xUnit*90,ylims(1,2)-yUnit*15,statsText,'HorizontalAlignment','Center','Interpreter','none','FontSize',8);
    end

    ylabel(x_label,'Interpreter','none');
    xlabel(y_label,'Interpreter','none');
else
    FullScatter(y,x,y_label,x_label,ScatterTitle,statsText) ; 
end
end
