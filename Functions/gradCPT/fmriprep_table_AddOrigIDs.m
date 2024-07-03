function [fmriprep_table] = fmriprep_table_AddOrigIDs(ExpDir,fmriprep_table)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ExpNames=unique(fmriprep_table.experiment);
for i = 1:length(ExpNames) 
    ExpInd=single(ismember(fmriprep_table.experiment,ExpNames{i,1}));
    subStruct=tdfread([ExpDir,ExpNames{i,1},'/participants.tsv'],'\t'); 
    BIDs_IDs=subStruct.participants_id;
    Sub_IDs=subStruct.sub_id;
    
    for j = 1:size(Sub_IDs,1)
        BIDs_ID = strrep(BIDs_IDs(j,:),' ','');
        BIDs_ID = strrep(BIDs_ID,'sub-','');
        Sub_ID = strrep(Sub_IDs(j,:),' ','');
        SubInd=(ExpInd + single(ismember(fmriprep_table.sub,BIDs_ID)))==2;
        fmriprep_table.orig_IDs(SubInd,1)={Sub_ID};
        fmriprep_table.TR(SubInd,1)=subStruct.TR(j,1);
    end
end    
end

