function [fmriprep_table] = fMRIPrep_Table_AddStartInds(fmriprep_table)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addVarNames={'sub_session_id','numRuns_bySub','numRuns_bySes'};
tableNames=fmriprep_table.Properties.VariableNames;
if any(ismember(addVarNames,tableNames))
    fmriprep_table(:,addVarNames)=[];
end
sub_session_id=strcat(fmriprep_table.sub,'_',num2str(fmriprep_table.session));
numRuns_bySes=~strcmpi([{''};sub_session_id(1:end-1,1)],sub_session_id);
b=find(numRuns_bySes);
c=[diff(b);length(numRuns_bySes)-b(end)+1];
numRuns_bySes=single(numRuns_bySes);
numRuns_bySes(b)=c;
fmriprep_table.numRuns_bySes=numRuns_bySes;
fmriprep_table.sub_session_id=sub_session_id;

subIDs=fmriprep_table.sub;
numRuns_bySub=~strcmpi([{''};subIDs(1:end-1,1)],subIDs);
b=find(numRuns_bySub);
c=[diff(b);length(numRuns_bySub)-b(end)+1];
numRuns_bySub=single(numRuns_bySub);
numRuns_bySub(b)=c;
fmriprep_table.numRuns_bySub=numRuns_bySub;
end

