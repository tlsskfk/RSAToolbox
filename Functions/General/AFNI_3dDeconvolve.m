 function [ runCode ] = AFNI_3dDeconvolve(inputDataPath,stimInfo,outFilePath,varargin )
%Written by David Rothlein
 %UNTITLED Summary of this function goes here
%   Detailed explanation goes here

runCode=[];
AM2_stimInfo = VariableSetter( 'AM2_stimInfo',[],varargin);
BaseLine_stimInfo = VariableSetter( 'BaseLine_stimInfo',[],varargin);
allzero_OK = VariableSetter( 'allzero_OK',1,varargin);
jobs = VariableSetter( 'jobs',1,varargin);
overwrite = VariableSetter( 'overwrite',1,varargin);
GOFORIT = VariableSetter( 'GOFORIT',4,varargin);
mask = VariableSetter( 'mask',[],varargin);
censor = VariableSetter( 'censor',[],varargin);
polort = VariableSetter( 'polort',[],varargin);
float = VariableSetter( 'float',1,varargin);
local_times = VariableSetter( 'local_times',1,varargin);
HDR = VariableSetter( 'HDR','GAM',varargin);
TR = VariableSetter( 'TR',[],varargin);
quiet = VariableSetter( 'quiet',1,varargin);
tout = VariableSetter( 'tout',1,varargin);
nofullf_atall = VariableSetter( 'nofullf_atall',1,varargin);
residOut = VariableSetter( 'residOut',[],varargin);
noFDR = VariableSetter( 'noFDR',1,varargin);
noData = VariableSetter( 'noData',0,varargin);
numVol = VariableSetter( 'numVol',[],varargin);
StartVol = VariableSetter( 'StartVol',[],varargin);
Use3DReml = VariableSetter( 'Use3DReml',0,varargin);
if iscell(TR)
    TR=TR{1,1};
end
labelcount=1;
filecount=1;
numEvents=size(stimInfo,1);
numAM2=size(AM2_stimInfo,1);
numBaseLine=size(BaseLine_stimInfo,1);
numStimt=numEvents+numAM2+numBaseLine;

%runCode=[runCode,'tcsh',newline,newline];
runCode=[runCode,'3dDeconvolve \',newline];

if allzero_OK==1
    runCode=[runCode,char(9),'-allzero_OK ',' \',newline];
end

if jobs>1
    runCode=[runCode,char(9),'-jobs ',num2str(jobs), ' \',newline];
end

if ~isempty(TR)
    runCode=[runCode,char(9),'-force_TR ',num2str(TR), ' \',newline];
end

if ~isempty(polort)
    runCode=[runCode,char(9),'-polort ',num2str(polort), ' \',newline];
end

if overwrite==1
    runCode=[runCode,char(9),'-overwrite ',' \',newline];
end

if GOFORIT>1
    runCode=[runCode,char(9),'-GOFORIT ',num2str(GOFORIT), ' \',newline];
end

if noData == 0
    runCode=[runCode,char(9),'-input ',inputDataPath, ' \',newline];
else
    runCode=[runCode,char(9),'-nodata ',num2str(numVol),' ',num2str(TR), ' \',newline];
    if ~isempty(StartVol)
        runCode=[runCode,char(9),'-concat ''1D: ',num2str(StartVol(:)'), ''' \',newline];
    end
end

if ~isempty(mask)
    runCode=[runCode,char(9),'-mask ',mask, ' \',newline];
end

if ~isempty(censor)
    runCode=[runCode,char(9),'-censor ',censor, ' \',newline];
end

if float==1
    runCode=[runCode,char(9),'-float ',' \',newline];
end

runCode=[runCode,char(9),'-num_stimts ',num2str(numStimt),' \',newline];

if local_times==1
    runCode=[runCode,char(9),'-local_times ',' \',newline];
end

if ~isempty(BaseLine_stimInfo)
    for i = 1:numBaseLine
        runCode=[runCode,char(9),char(9),'-stim_label ',num2str(labelcount),' "',BaseLine_stimInfo{i,1},'" \',newline];
        labelcount=labelcount+1;  
    end
    if ~isempty(stimInfo)
        for i = 1:numBaseLine
            runCode=[runCode,char(9),char(9),'-stim_file ',num2str(filecount),' "',BaseLine_stimInfo{i,2},'" ','-stim_base ',num2str(filecount),' \',newline];
            filecount=filecount+1;
        end
    else
        for i = 1:numBaseLine
            runCode=[runCode,char(9),char(9),'-stim_file ',num2str(filecount),' "',BaseLine_stimInfo{i,2},'" ',' \',newline];
            filecount=filecount+1;
        end
    end
end
if ~isempty(stimInfo)
    for i = 1:numEvents
        runCode=[runCode,char(9),char(9),'-stim_label ',num2str(labelcount),' "',stimInfo{i,1},'" \',newline];
        labelcount=labelcount+1;
    end

    for i = 1:numEvents
        runCode=[runCode,char(9),char(9),'-stim_times ',num2str(filecount),' "',stimInfo{i,2},'" ',HDR,' \',newline];
        filecount=filecount+1;
    end
end

if ~isempty(AM2_stimInfo)
    for i = 1:numAM2
        runCode=[runCode,char(9),char(9),'-stim_label ',num2str(labelcount),' "',AM2_stimInfo{i,1},'" \',newline];
        labelcount=labelcount+1;
    end
end
if ~isempty(AM2_stimInfo)
    for i = 1:numAM2
        runCode=[runCode,char(9),char(9),'-stim_times_AM2 ',num2str(filecount),' "',AM2_stimInfo{i,2},'" ',HDR,' \',newline];
        filecount=filecount+1; 
    end
end
if tout==1
    runCode=[runCode,char(9),'-tout ',' \',newline];
end

if nofullf_atall==1
    runCode=[runCode,char(9),'-nofullf_atall ',' \',newline];
end

if noFDR==1
    runCode=[runCode,char(9),'-noFDR ',' \',newline];
end

if quiet==1
    runCode=[runCode,char(9),'-quiet ',' \',newline];
end

if ~isempty(residOut)
    runCode=[runCode,char(9),'-errts ',residOut,' \',newline];
end
if Use3DReml==1
    runCode=[runCode,char(9),'-x1D_stop ',' \',newline];
end 
runCode=[runCode,char(9),'-bucket ',outFilePath,newline];

    
if Use3DReml == 1
    runCode=[newline,runCode,'3dTcat -prefix ',[outFilePath,'AllRuns '],strrep(inputDataPath,'.HEAD',''),newline];
    runCode=[newline,runCode,'3dREMLfit  \',newline];
    runCode=[runCode,char(9),'-input ',[outFilePath,'AllRuns+tlrc.HEAD'], ' \',newline];
    runCode=[runCode,char(9),'-matrix ',outFilePath, '.xmat.1D \',newline];
    if GOFORIT>1
        runCode=[runCode,char(9),'-GOFORIT ',num2str(GOFORIT), ' \',newline];
    end
    if tout==1
        runCode=[runCode,char(9),'-tout ',' \',newline];
    end    
    if noFDR==1
        runCode=[runCode,char(9),'-noFDR ',' \',newline];
    end  
    if overwrite==1
        runCode=[runCode,char(9),'-overwrite ',' \',newline];
    end    
    runCode=[runCode,char(9),'-RBuck ',outFilePath,newline];
end
    