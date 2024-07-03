function [FileNames] = AFNIonsetMaker_stimfile(TimeCourses,SaveDir,SaveName,AltDir)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin==3
        AltDir=SaveDir;
    end

    if ~iscell(TimeCourses)
        TimeCourses={TimeCourses};
    end
    TimeCourses=TimeCourses(:);
    numRuns=length(TimeCourses);

    AllData=[];
    for i = 1:numRuns
        AllData=[AllData;TimeCourses{i,1}(:)];          
    end

    FileName=[SaveDir,'/',SaveName,'.txt'];
    FileName=strrep(FileName,'//','/');
    AFNI_Onset_Maker_stimfile( FileName,AllData);
    FileNames=[AltDir,'/',SaveName,'.txt']; 
    FileNames=strrep(FileNames,'//','/');
end

function AFNI_Onset_Maker_stimfile( SaveName,Data)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    fileID = fopen(SaveName,'w');
    fprintf(fileID,'%.6f\r\n',Data);
    fclose(fileID);
end