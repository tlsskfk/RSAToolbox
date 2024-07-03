function [Results] = MetaGroupvarSummary(ExperimentsDir)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    addData=1;
    addCount=1;
    
    ExperimentsDir=strrep(ExperimentsDir,'\','/');
    GroupDir=strrep(ExperimentsDir,'Experiments/','GroupSummaries/');
    GroupDir=strrep(GroupDir,'//','/'); 
    LoadSummary=uiNameSelect({'Yes','No'},'Load existing data summary?',1);
    if strcmpi(LoadSummary,'Yes')
        [~,SummaryFileNames]=getFolderAndFileNames(GroupDir);
        SummaryFileNames=uiNameSelect(SummaryFileNames,'Select data summary file',1); 
        if iscell(SummaryFileNames)
            SummaryFileNames=SummaryFileNames{1,1};
        end   
        load([GroupDir,SummaryFileNames],'Results','ParcelNames'); 
        addCount=2;
    else
        SummaryFileNames='';
    end
    if ~exist(GroupDir,'file')
        mkdir(GroupDir);
    end   
        
    [SaveName] = uiEnterName(SummaryFileNames,'Enter save name for data summary file');
        
    while addData==1
        [tempResults,NewParcelNames] = GroupVarSummary(ExperimentsDir,'NoSave',1); 
        [UsePrefix] = uiEnterName(['Group',num2str(addCount),'_'],'Enter group prefix');
        tempResults=TableCellAppendName(tempResults,UsePrefix,'');
        addData=uiNameSelect({'Yes','No'},'Add more data?',1);
        if addCount == 1
            Results=tempResults;
            ParcelNames=NewParcelNames;
        else
            [Results,ParcelNames]=TableCellCat(Results,tempResults,ParcelNames,NewParcelNames);
        end
        if strcmpi(addData,'Yes')
            addData=1;
            addCount=addCount+1;
        else
            addData=0;
        end
        save([GroupDir,SaveName],'Results','ParcelNames');
    end
end

function [CatCell,ParcelNames]=TableCellCat(OldCell,NewCell,OldParcelNames,NewParcelNames)
    OldParcelNames=OldParcelNames(:);
    NewParcelNames=NewParcelNames(:);
    ParcelNames=OldParcelNames;
    CatCell=OldCell;
    for i = 1:length(NewParcelNames)
        numCol=size(NewCell,2);
        parcelName=NewParcelNames{i,1};
        if any(find(ismember(OldParcelNames,parcelName)))
            rowInd=find(ismember(OldParcelNames,parcelName));
            for j = 1:numCol
                CatCell{rowInd,j}=[OldCell{rowInd,j},NewCell{i,j}];
            end
        else
            oldCol=size(CatCell,2);
            if numCol > oldCol
                CatCell=[CatCell,cell(size(CatCell,1),numCol-oldCol)];
            elseif numCol < oldCol
                NewCell=[NewCell,cell(size(NewCell,1),oldCol-numCol)];
            end
            CatCell=[CatCell;NewCell(i,:)];
            ParcelNames=[ParcelNames;{parcelName}];
        end
    end
end

function TableCell=TableCellAppendName(TableCell,Prefix,Suffix)
    numRow=size(TableCell,1);
    numCol=size(TableCell,2);
    for i = 1:numRow
        for j = 1:numCol
            if ~isempty(TableCell{i,j})
                VarNames=TableCell{i,j}.Properties.VariableNames(:);
                for k = 1:length(VarNames)
                    VarNames{k,1}=[Prefix,VarNames{k,1},Suffix];
                end
                TableCell{i,j}.Properties.VariableNames=VarNames;
            end
        end
    end
    
end