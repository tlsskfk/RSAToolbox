function [b] = FullScatter(y,x,y_label,x_label,ScatterTitle,statsText)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    y_label='';
    x_label='';
    ScatterTitle='';
    statsText='';
end
if nargin==3
    x_label='';    
    statsText='';
    ScatterTitle='';    
end
if nargin==4  
    statsText='';
    ScatterTitle='';    
end
if nargin==5  
    statsText='';   
end
figure;
scatter(y,x)
hold on
title(ScatterTitle,'Interpreter','none');
b=regress(x,[y,ones(length(y),1)]);
slope=b(1,1);
Int=b(2,1);
xRange=[min(y),max(y)];
yRange=[Int+xRange(1,1)*slope,Int+xRange(1,2)*slope];
plot(xRange,yRange,'k','LineWidth',2);
text(mean(xRange(:)),min(x),statsText,'HorizontalAlignment','Center','Interpreter','none')
ylabel(x_label,'Interpreter','none');
xlabel(y_label,'Interpreter','none');
end

