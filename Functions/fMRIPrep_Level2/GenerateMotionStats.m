function [fmriprep_table] = GenerateMotionStats(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%Template function for data processing from the fmriprep_table

%% Set default values using variable setter function
%Subject or run level analysis. Will prompt request.
[SubjectOrRun] = VariableSetter('SubjectOrRun','Subject',varargin);
[TR] = VariableSetter('TR',[],varargin);
[VaryingTR] = VariableSetter('VaryingTR',[],varargin);
%% Compute or set initial parameters and names for the analysis
%Compute total number of runs and number of participants 
TotalRuns=height(fmriprep_table);
dataInd=find(fmriprep_table.run==1)';
numSS=length(dataInd);

%Set analysis type and analysis name. These values will be used when saving
%files and adding colunms to the BIDs table.

%Select whether to run analysis for each run seperately or combine data 
%across runs and perform analysis by SS. Treats sessions seperately (but will
%include an option to combine across session in future).
if isempty(SubjectOrRun)
    SingleSelect=1; %Allows only a single value to be selected.
    [SubjectOrRun] = uiNameSelect({'Subject','Run'},'Perform analysis by subject or by run:',SingleSelect);
end
% if strcmpi(SubjectOrRun,'Subject')
%     bySS=1;
%     useIndicies=dataInd;
% else
%     bySS=0;
%     useIndicies=[1:TotalRuns];
% end
if strcmpi(SubjectOrRun,'Subject')
    bySS=1;
    if length(unique(fmriprep_table.session))>1
        useIndicies=find(fmriprep_table.numRuns_bySes)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySes;        
    else
        useIndicies=find(fmriprep_table.numRuns_bySub)';
        fmriprep_table.numRuns=fmriprep_table.numRuns_bySub;
    end
else
    bySS=0;
    useIndicies=[1:TotalRuns];
end
if ~any(ismember(fmriprep_table.Properties.VariableNames,'TR'))
    InputTR=1;
else
    InputTR=0;
end
if isempty(TR) && InputTR==1
    [VaryingTR] = uiNameSelect({'Yes','No'},'Does TR vary across participant?:',1);
    if strcmpi(VaryingTR,'No')
        TR=uiEnterName('','Enter TR (sec):');
        TR=str2num(TR);
        VaryingTR=0;
    else
        VaryingTR=1;
    end
else    
    VaryingTR=0;
end
%% BIDsTable loop: Iterate through fmriprep_table and perform analysis
iniPercentComplete=0; %Used to display progress

for dataInd=useIndicies
    tic
    %% Display progress
    PercentComplete=round((dataInd/TotalRuns)*100);
    if PercentComplete>iniPercentComplete
        disp([num2str(PercentComplete),'% Complete. Generate Motion Stats']);
        iniPercentComplete=PercentComplete;
    end  
    
    %% Initialize input data for loading
    if bySS==1
        numRuns=fmriprep_table.numRuns(dataInd,1);
    else
        numRuns=1;
    end
    count=1;
    Rot_bySS=[];
    Trans_bySS=[];
    nsso_bySS=0;
    fd_bySS=[];
    %% load input data 
    % If by Subject, iterate through runs and place data in cell
    for run=1:numRuns
        loadInd=dataInd+run-1;   
        %% set load paths and variable names
        LoadVars={'confounds_regressors'};
        try
            LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.funcDir{loadInd,1},fmriprep_table.(LoadVars{1,1}){loadInd,1}];
        catch
            LoadVars={'confounds_timeseries'};
            LoadPaths{1,1}=[ExperimentsDir,fmriprep_table.funcDir{loadInd,1},fmriprep_table.(LoadVars{1,1}){loadInd,1}];
        end
        try
            [ ConfoundTCs ] = ConfoundLoad( LoadPaths{1,1},1 ); 
            fmriprep_table.numVol{loadInd,1}=size(ConfoundTCs,1);
            if InputTR==1
                if VaryingTR==1
                    subID=[fmriprep_table.sub{loadInd,1},'_run',num2str(fmriprep_table.run(loadInd,1)),'_',fmriprep_table.task{loadInd,1}];
                    TR=uiEnterName('',['Enter TR (sec) for: ',newline,subID]);
                    TR=str2num(TR);
                end
                fmriprep_table.TR{loadInd,1}=TR;
            end
            try
                fmriprep_table.fmriDur{loadInd,1}=fmriprep_table.numVol{loadInd,1}*fmriprep_table.TR{loadInd,1};
            catch
                fmriprep_table.fmriDur{loadInd,1}=fmriprep_table.numVol{loadInd,1}*fmriprep_table.TR(loadInd,1);
            end
            rot=[];
            trans=[];
            nsso=0;
            fd=[];
            for j = 1:width(ConfoundTCs)
                if contains(ConfoundTCs(:,j).Properties.VariableNames,'rot')
                    rot=[rot;table2array(ConfoundTCs(:,j))];
                elseif contains(ConfoundTCs(:,j).Properties.VariableNames,'trans')
                    trans=[trans;table2array(ConfoundTCs(:,j))];
                elseif contains(ConfoundTCs(:,j).Properties.VariableNames,'non_steady_state')
                    nsso=nansum(table2array(ConfoundTCs(:,j)),1);
                elseif contains(ConfoundTCs(:,j).Properties.VariableNames,'framewise_displacement')
                    fd=table2array(ConfoundTCs(:,j));
                end
            end
            fmriprep_table.MotionByRun_rot_mean(loadInd,1)=nanmean(rot(:));
            fmriprep_table.MotionByRun_rot_max(loadInd,1)=nanmax(rot(:));
            fmriprep_table.MotionByRun_trans_mean(loadInd,1)=nanmean(trans(:));
            fmriprep_table.MotionByRun_trans_max(loadInd,1)=nanmax(trans(:)); 
            fmriprep_table.MotionByRun_nsso_sum(loadInd,1)=nsso;
            fmriprep_table.MotionByRun_fd_mean(loadInd,1)=nanmean(fd(:));
            fmriprep_table.MotionByRun_fd_max(loadInd,1)=nanmax(fd(:));            
            Rot_bySS=[Rot_bySS;rot];
            Trans_bySS=[Trans_bySS;trans];
            nsso_bySS=nsso_bySS+nsso;
            fd_bySS=[fd_bySS;fd];
            count=count+1;
        catch
            disp(['Skipping run-- input file or variable doesnt exist: ',LoadVars{1,1},newline,LoadPaths{1,1}]);
            continue
        end          
    end      
    if count==1
        disp(['Skipping subject-- no input files or variables exist: ',LoadVars{1,1},newline,LoadPaths{1,1}]);
        continue
    end  
    
    %% correct for non steady state outliers
    fmriprep_table.MotionBySS_fd_mean(dataInd,1)=nanmean(fd_bySS(:));
    fmriprep_table.MotionBySS_fd_max(dataInd,1)=nanmax(fd_bySS(:));
    fmriprep_table.MotionBySS_rot_mean(dataInd,1)=nanmean(Rot_bySS(:));
    fmriprep_table.MotionBySS_rot_max(dataInd,1)=nanmax(Rot_bySS(:));    
    fmriprep_table.MotionBySS_trans_mean(dataInd,1)=nanmean(Trans_bySS(:));
    fmriprep_table.MotionBySS_trans_max(dataInd,1)=nanmax(Trans_bySS(:)); 
    fmriprep_table.MotionBySS_nsso_sum(dataInd,1)=nsso_bySS;       
    toc
end
end 

function [ ConfoundRegressors ] = ConfoundLoad( FileName,combineNSSO )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ConfoundStruct=tdfread(FileName,'\t');
VarNames=fieldnames(ConfoundStruct);
for i=1:length(VarNames)
    VarName=VarNames{i,1};
    if ~isnumeric(ConfoundStruct.(VarName))
        newMat=nan(size(ConfoundStruct.(VarName),1),1);
        removeInd=ismember(ConfoundStruct.(VarName)(:,1),'n');
        ConfoundStruct.(VarName)(removeInd,:)=[];
        tempMat=str2num(ConfoundStruct.(VarName));
        newMat(removeInd==0)=tempMat;
        ConfoundStruct.(VarName)=newMat;
    end
    if combineNSSO==1   
        if contains(VarName,'non_steady_state_outlier')
            if isfield(ConfoundStruct,'non_steady_state_outlier')
                ConfoundStruct.non_steady_state_outlier=ConfoundStruct.non_steady_state_outlier+ConfoundStruct.(VarName);
            else
                ConfoundStruct.non_steady_state_outlier=ConfoundStruct.(VarName);
            end
            ConfoundStruct=rmfield(ConfoundStruct,VarName);
        end
    end
end
ConfoundStruct.Linear=[1:length(ConfoundStruct.(VarNames{1,1}))]'/length(ConfoundStruct.(VarNames{1,1}));
ConfoundRegressors=struct2table(ConfoundStruct);

end