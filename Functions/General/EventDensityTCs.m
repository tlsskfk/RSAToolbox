function [ EventDensity,preEvent,postEvent ] = EventDensityTCs( EventTrials,gauss_mu,gauss_sigma)
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
EventTrials=EventTrials(:);
nanInd=isnan(EventTrials);
EventTrials(nanInd)=0;
numEvents=sum(EventTrials(:));
numTrials=length(EventTrials);
if numEvents > 0   
    gaussWindow = gauss_sigma * 4;
    Gauss = normpdf([0:gaussWindow],gauss_mu,gauss_sigma);
    postEvent=conv(EventTrials,Gauss);
    preEvent=flip(conv(flip(EventTrials,1),Gauss),1);
    preEvent=[preEvent(gaussWindow+1:end);zeros(gaussWindow,1)];
    preEvent=preEvent(1:numTrials);
    postEvent=postEvent(1:numTrials);
    EventDensity=postEvent+preEvent;
    EventDensity(EventTrials==1)=nan;
    EventDensity=fillmissing(EventDensity,'linear');    
    preEvent=scale01(preEvent);
    postEvent=scale01(postEvent);
    EventDensity=scale01(EventDensity);
    preEvent(nanInd)=nan;
    postEvent(nanInd)=nan;
    EventDensity(nanInd)=nan;       
else
    EventDensity=zeros(length(EventTrials),1);
    preEvent=zeros(length(EventTrials),1);
    postEvent=zeros(length(EventTrials),1);    
end

end

