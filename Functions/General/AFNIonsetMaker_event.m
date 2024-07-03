function [FileNames] = AFNIonsetMaker_event(Events,OnsetTimes,SaveDir,SaveName,AltDir)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin==4
    AltDir=SaveDir;
end
    
if ~iscell(Events)
    Events={Events};
end
Events=Events(:);
numRuns=length(Events);

if ~iscell(OnsetTimes)
    OnsetTimes={OnsetTimes};
end
OnsetTimes=OnsetTimes(:);

for i = 1:numRuns
    Event=Events{i,1};
    OnsetTime=OnsetTimes{i,1};
    EventTimes=OnsetTime(Event==1,:)';
    if i == 1
        NewLine=0;
        FileName=[SaveDir,'/',SaveName,'.txt'];
        FileName=strrep(FileName,'//','/');
        AFNI_Onset_Maker( FileName,EventTimes,NewLine);
        FileNames=[AltDir,'/',SaveName,'.txt']; 
        FileNames=strrep(FileNames,'//','/');
    else
        NewLine=1;
        FileName=[SaveDir,'/',SaveName,'.txt'];
        FileName=strrep(FileName,'//','/');
        AFNI_Onset_Maker( FileName,EventTimes,NewLine);
        FileNames=[AltDir,'/',SaveName,'.txt']; 
        FileNames=strrep(FileNames,'//','/');
    end                
end

