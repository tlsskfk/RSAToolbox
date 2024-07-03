
function [OutputTable,SaveName,AnalysisParameters] = GroupTableMaker(ExperimentsDir,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[TableNames] = VariableSetter('TableNames',[],varargin);
[AnalysisTypes] = VariableSetter('AnalysisTypes',[],varargin);
[AnalysisNames] = VariableSetter('AnalysisNames',[],varargin);
[ParcellationNames] = VariableSetter('ParcellationNames',[],varargin);
[LoadVariableName] = VariableSetter('LoadVariableName',[],varargin);
[LoadFileName] = VariableSetter('LoadFileName',[],varargin);
[TableRowNames] = VariableSetter('TableRowNames',[],varargin);
[TableVariableNames] = VariableSetter('TableVariableNames',[],varargin);
[tempTableVarNames] = VariableSetter('tempTableVarNames',[],varargin);
[UseVarNames] = VariableSetter('UseVarNames',[],varargin);
[NewVarOrder] = VariableSetter('NewVarOrder',[],varargin);
[NewVarNames] = VariableSetter('NewVarNames',[],varargin);
AnalysisParameters = VariableSetter('AnalysisParameters',[],varargin);

if isempty(AnalysisParameters)
    SingleSelect=1; %Allows only a single value to be selected.
    [LoadParams] = uiNameSelect({'Yes','No'},'Load existing parameters?',SingleSelect);   
    if strcmpi(LoadParams,'Yes')
        try
            ExperimentsDir=strrep(ExperimentsDir,'\','/');
            GroupTableDir=strrep(ExperimentsDir,'Experiments/','GroupTables/');
            GroupTableDir=strrep(GroupTableDir,'//','/'); 
            TableTypeName=getFolderAndFileNames(GroupTableDir);
            TableTypeName=uiNameSelect(TableTypeName,'Select tabletype',1); 
            if iscell(TableTypeName)
                TableTypeName=TableTypeName{1,1};
            end
            GroupTableDir=[GroupTableDir,TableTypeName,'/'];
            [~,TableTypeFileNames]=getFolderAndFileNames(GroupTableDir);
            TableTypeFileNames=uiNameSelect(TableTypeFileNames,'Select table file',1); 
            if iscell(TableTypeFileNames)
                TableTypeFileNames=TableTypeFileNames{1,1};
            end   
            TableTypeAnalysisParams=load([GroupTableDir,TableTypeFileNames],'AnalysisParameters');
            AnalysisParameters=TableTypeAnalysisParams.AnalysisParameters;
        catch
            disp('No Existing Parameters!')
        end
    end
end

if ~isempty(AnalysisParameters)
    ParamNames=fieldnames(AnalysisParameters);
%     for i = 1:length(ParamNames)
%         if isempty(AnalysisParameters.(ParamNames{i,1}))
%             AnalysisParameters.(ParamNames{i,1})='NA_EMPTY';
%         end
%     end    
    SingleSelect=0; %Allows only a single value to be selected.
    [EditParams] = uiNameSelect(fieldnames(AnalysisParameters),'Select analysis parameters to edit:',SingleSelect);    
    if ~isempty(EditParams)
        for i = 1:length(EditParams)
            AnalysisParameters.(EditParams{i,1})=[];
        end
    end
    TableNames=AnalysisParameters.TableNames;
    AnalysisTypes=AnalysisParameters.AnalysisTypes;
    AnalysisNames=AnalysisParameters.AnalysisNames;
    ParcellationNames=AnalysisParameters.ParcellationNames;
    LoadVariableName=AnalysisParameters.LoadVariableName;
    LoadFileName=AnalysisParameters.LoadFileName;
    TableRowNames=AnalysisParameters.TableRowNames;
    TableVariableNames=AnalysisParameters.TableVariableNames;
    tempTableVarNames=AnalysisParameters.tempTableVarNames;
    UseVarNames=AnalysisParameters.UseVarNames;
    NewVarOrder=AnalysisParameters.NewVarOrder;
    NewVarNames=AnalysisParameters.NewVarNames; 
else
    AnalysisParameters=struct;
end


%% Compile files to use to make table
[AllPaths_Parcel,AllPaths_NoParcel,FolderNames_Parcel,FolderNames_NoParcel,GroupDir,AnalysisNames,AnalysisTypes,TableNames,LoadFileName,ParcellationNames] = GroupDirSearch(ExperimentsDir,'TableNames',TableNames,'AnalysisTypes',AnalysisTypes,'AnalysisNames',AnalysisNames,'ParcellationNames',ParcellationNames,'FileNames',LoadFileName,'GroupType','GroupAnalysis');
if isempty(AllPaths_Parcel)
    UsePaths=AllPaths_NoParcel;
    UseFolderNames=FolderNames_NoParcel;
    Parcel=0;
else
    UsePaths=AllPaths_Parcel;
    UseFolderNames=FolderNames_Parcel;
    Parcel=1;    
end
SaveDir=strrep(GroupDir,'GroupAnalysis','GroupTables');

UseIDVarNames=[];
FigureFolderName=[];
FigureFileName=[];
for i = 1:size(UseFolderNames,2)
    if height(unique(UseFolderNames(:,i)))>1
        UseIDVarNames=[UseIDVarNames,UseFolderNames.Properties.VariableNames(i)];
        FigureFolderName=[FigureFolderName,UseFolderNames.Properties.VariableNames{i},'_'];
    else
        FigureFileName=[FigureFileName,UseFolderNames.(UseFolderNames.Properties.VariableNames{i}){1,1},'_'];
    end
end
FigureFolderName(1,end)='/';
FigureFolderName=[SaveDir,FigureFolderName];
FigureFileName=strrep(FigureFileName,'.mat','');
if ~exist(FigureFolderName,'file')
    mkdir(FigureFolderName);
end
if ~iscell(UseIDVarNames)
    UseIDVarNames={UseIDVarNames};
end
IDFolderNames=UseFolderNames(:,UseIDVarNames);
if isempty(UseVarNames)
    UseVarNames=join(table2cell(IDFolderNames),2);
    UseVarNames=strrepCell(UseVarNames,' ','_');
    UseVarNames=strrepCell(UseVarNames,'.mat','');
end
for i = 1:length(UseVarNames)
   if length(UseVarNames{i,1}) > 63
       UseVarNames{i,1}=UseVarNames{i,1}(1,end-62:end);
   end
end    
%% Select which variable to load from mat file
if isempty(LoadVariableName)
    for i = 1:length(UsePaths)
        variableInfo = who('-file', AllPaths_Parcel{i,1});
        LoadVariableName=unique([LoadVariableName;variableInfo]);
    end
    [LoadVariableName] = uiNameSelect(LoadVariableName,'Select variable to load:' ,1);
end
if iscell(LoadVariableName)
    LoadVariableName=LoadVariableName{1,1};
end
FigureFileName=[FigureFileName,LoadVariableName,'_'];
%% Compile Table Variable Names and Row Names
SkipSelect=0;
if any(isempty(TableRowNames)) || any(strcmpi(TableRowNames,'All'))
    for i = 1:length(UsePaths)
        UseTable=load(UsePaths{i,1},LoadVariableName);
        UseTable=UseTable.(LoadVariableName);
        if isempty(UseTable.Properties.RowNames)
            SkipSelect=1;
            break
        elseif i==1 && strcmpi(TableRowNames,'All')
            TableRowNames=[];
            SkipSelect=1;
        end           
        TableRowNames=unique([TableRowNames;UseTable.Properties.RowNames]);
    end 
    if SkipSelect==0
        [TableRowNames] = uiNameSelect(TableRowNames,'Select rownames to include:' ,0);
    end
end
tempTableVariableNames=[];
if isempty(TableRowNames) && isempty(tempTableVarNames)
    for i = 1:length(UsePaths)
        UseTable=load(UsePaths{i,1},LoadVariableName);
        UseTable=UseTable.(LoadVariableName);            
        tempTableVariableNames=unique([tempTableVariableNames;UseTable.Properties.VariableNames(:)]);
    end  
    [tempTableVarNames] = uiNameSelect(tempTableVariableNames,'Select variable to set as row names:' ,1);
end
SkipSelect=0;
if isempty(TableRowNames) || any(strcmpi(TableRowNames,'All'))
    for i = 1:length(UsePaths)
        UseTable=load(UsePaths{i,1},LoadVariableName);
        UseTable=UseTable.(LoadVariableName);  
        UseTable.Properties.RowNames=UseTable.(tempTableVarNames);
        if i==1 && strcmpi(TableRowNames,'All')
            TableRowNames=[];
            SkipSelect=1;
        end                  
        TableRowNames=unique([TableRowNames;UseTable.Properties.RowNames]);
    end 
    if SkipSelect==0
        [TableRowNames] = uiNameSelect(TableRowNames,'Select rownames to include:' ,0);
    end
end 

if isempty(TableVariableNames)  
    for i = 1:length(UsePaths)
        UseTable=load(UsePaths{i,1},LoadVariableName);
        UseTable=UseTable.(LoadVariableName);            
        TableVariableNames=unique([TableVariableNames;UseTable.Properties.VariableNames(:)]);
    end 
    [TableVariableNames] = uiNameSelect(TableVariableNames,'Select variables to include:' ,0);
elseif ~iscell(TableVariableNames)
    TableVariableNames={TableVariableNames};
end

if isempty(NewVarOrder)
    [NewVarOrder,NewVarNames] = uiTableVarNames(UseVarNames)   ; 
end

OutputTable=[];
for i = 1:length(UsePaths)
    UseTable=load(UsePaths{i,1},LoadVariableName);
    UseTable=UseTable.(LoadVariableName);   
    if ~isempty(tempTableVarNames)
        UseTable.Properties.RowNames=UseTable.(tempTableVarNames);
    end
    tempTable=UseTable(TableRowNames,TableVariableNames);
    if length(TableVariableNames)==1
        tempTable.Properties.VariableNames=UseVarNames(i,1);
        OutputTable=[OutputTable,tempTable];
    else
        tempTable=rows2vars(tempTable);
        tempTable.Properties.RowNames=tempTable.OriginalVariableNames;
        tempTable.OriginalVariableNames=[];
        tempTable.Properties.VariableNames=UseVarNames(i,1);
        OutputTable=[OutputTable,tempTable];
    end
end  
temVarNames=OutputTable.Properties.VariableNames;
SelectInd=ismember(NewVarOrder,temVarNames);
NewVarOrder=NewVarOrder(1,SelectInd==1);
NewVarNames=NewVarNames(1,SelectInd==1);
OutputTable=OutputTable(:,NewVarOrder);
OutputTable.Properties.VariableNames=NewVarNames;
if length(TableVariableNames)==1
    FigureFileName=[FigureFileName,TableVariableNames,'_',genDateString];
else
    FigureFileName=[FigureFileName,TableRowNames,'_',genDateString];
end

AnalysisParameters.TableNames=TableNames;
AnalysisParameters.AnalysisTypes=AnalysisTypes;
AnalysisParameters.AnalysisNames=AnalysisNames;
AnalysisParameters.ParcellationNames=ParcellationNames;
AnalysisParameters.LoadVariableName=LoadVariableName;
AnalysisParameters.LoadFileName=LoadFileName;
AnalysisParameters.TableRowNames=TableRowNames;
AnalysisParameters.TableVariableNames=TableVariableNames;
AnalysisParameters.tempTableVarNames=tempTableVarNames;
AnalysisParameters.UseVarNames=UseVarNames;
AnalysisParameters.NewVarOrder=NewVarOrder;
AnalysisParameters.NewVarNames=NewVarNames; 

SaveName=strrep(strjoin([FigureFolderName,FigureFileName]),' ','');
save(SaveName,'OutputTable','AnalysisParameters');

end
 
