function [FileNames] = AFNIonsetMaker_AM2TC(TCs,AM2Event,OnsetTimes,SaveDir,SaveName,AltDir)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin==5
    AltDir=SaveDir;
end
    
if ~iscell(TCs)
    TCs={TCs};
end
TCs=TCs(:);
numRuns=length(TCs);

if ~iscell(AM2Event)
    AM2Event={AM2Event};
end
AM2Event=AM2Event(:);

if ~iscell(OnsetTimes)
    OnsetTimes={OnsetTimes};
end
OnsetTimes=OnsetTimes(:);

for i = 1:numRuns
    AllData=[OnsetTimes{i,1},TCs{i,1}];
    Event=AM2Event{i,1};
    EventTimes=AllData(Event==1,:)';   
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

