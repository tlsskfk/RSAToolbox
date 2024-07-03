function [AllDone] = DeleteAnalysisNames(ExperimentsDir,fmriprep_table)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
AllDone=0;
[~,~,~,~,~,~,~,~,AllAnalysisTypes]=BIDsDirSearch(ExperimentsDir,fmriprep_table,'SubjectOrRun','Run','GetAnalysisTypes',1);
UseAnalysisTypes=uiNameSelect(AllAnalysisTypes,'Select analysis types to search',0);
if ~iscell(UseAnalysisTypes)
    UseAnalysisTypes={UseAnalysisTypes};
end
UseAnalysisNames=cell(1,length(UseAnalysisTypes));

for i =1:length(UseAnalysisTypes)
    try
        [~,~,~,~,~,~,~,AllAnalysisNames,~]=BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',UseAnalysisTypes{i},'SubjectOrRun','Run','GetAnalysisNames',1);
        UseAnalysisNames{1,i}=uiNameSelect(AllAnalysisNames,[UseAnalysisTypes{i},': Select analysis names to delete'],0);
    catch
        UseAnalysisNames{1,i}={[]};
    end
end

for i =1:length(UseAnalysisTypes)
    UseAnalysisType=UseAnalysisTypes{i};
    tempAnalysisNames=UseAnalysisNames{1,i};
    if ~iscell(tempAnalysisNames)
        tempAnalysisNames={tempAnalysisNames};
    end 
    try
        if isempty(tempAnalysisNames{1})
            continue
        end   
    catch
        continue
    end
    for j = 1:length(tempAnalysisNames)
        tempAnalysisName=tempAnalysisNames{j};
        try
            BIDsDirSearch(ExperimentsDir,fmriprep_table,'AnalysisType',UseAnalysisType,'AnalysisName',tempAnalysisName,'SubjectOrRun','Run','DeleteAnalysisName',1);
            disp(['Deleted - ',UseAnalysisType,' ',tempAnalysisName]);
        catch
            disp(['Error - ',UseAnalysisType,' ',tempAnalysisName]);
        end
    end
end
AllDone=1;
end