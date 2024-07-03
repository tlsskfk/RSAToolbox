function [LoadVarName,DataType,LoadVarFormat,CompileType,CompileFormat,ndim] = GetFileVariableNames(filePaths,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[LoadVarName] = VariableSetter('LoadVarName',[],varargin);
[LoadVarFormat] = VariableSetter('LoadVarFormat',[],varargin);
[DataType] = VariableSetter('DataType',[],varargin);
[CompileType] = VariableSetter('CompileType',[],varargin); %TableOrArray
[CompileFormat] = VariableSetter('CompileFormat',[],varargin); %VertOrMat
[ndim] = VariableSetter('ndim',[],varargin);
ParcelNames=filePaths.Properties.VariableNames(:);


if isempty(LoadVarName)
    for i = 1:size(ParcelNames,1)
        ParcelName=ParcelNames{i,1};
        for j = 1:size(filePaths.(ParcelName),1)
            try
                variableInfo = who('-file', filePaths.(ParcelName){j,1});
            catch
                continue
            end        
            LoadVarName=unique([LoadVarName;variableInfo]);
            break
        end
    end
    [LoadVarName] = uiNameSelect(LoadVarName,'Which variable do you want to compile?',1);
end


if length(ParcelNames) > 1
    DataType='ByParcellation';
elseif isempty(DataType)
    SingleSelect=1; %Allows only a single value to be selected.
    [DataType] = uiNameSelect({'ByParcellation','ByCoordinate','Other'},'How is the data organized?',SingleSelect);
end  

if isempty(LoadVarFormat) || isempty(CompileType) || isempty(CompileFormat) || isempty(ndim)
    for i = 1:size(ParcelNames,1)
        ParcelName=ParcelNames{i,1};
        tempData=[];
        for j = 1:size(filePaths.(ParcelName),1)
            try
                tempData=load(filePaths.(ParcelName){j,1},LoadVarName);
                tempData=tempData.(LoadVarName);
            catch
                continue
            end        
            break
        end
        if ~isempty(tempData)
            break
        end
    end
    if istable(tempData)
        LoadVarFormat='Table';
        if any(ismember(size(tempData),1))
            ndim=1;  
            if isempty(CompileType)
                [CompileType] = uiNameSelect({'Table','Array'},'Compile in table or array form?',1);
            end       
            if strcmpi(CompileType,'Array') && strcmpi(DataType,'ByParcellation') && isempty(CompileFormat)
                [CompileFormat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',1);
            else
                CompileFormat='Matrix';
            end
        else
            ndim=2;
            CompileType='Array';
            CompileFormat='Matrix';
        end
    elseif iscell(tempData)
        LoadVarFormat='Cell';
        CompileType='Array';
        if ndims(tempData)>2
            ndim=ndims(tempData);
            CompileFormat='Matrix';
        elseif any(ismember(size(tempData),1))
            if strcmpi(CompileType,'Array') && strcmpi(DataType,'ByParcellation') && isempty(CompileFormat)
                [CompileFormat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',1);
            else
                CompileFormat='Matrix';
            end
            ndim=1;
        else
            ndim=2;
            CompileFormat='Matrix';
        end
    else
        LoadVarFormat='Array';
        CompileType='Array';
        if ndims(tempData)>2
            ndim=ndims(tempData);
            CompileFormat='Matrix';
        elseif any(ismember(size(tempData),1))
            if strcmpi(CompileType,'Array') && strcmpi(DataType,'ByParcellation') && isempty(CompileFormat)
                [CompileFormat] = uiNameSelect({'Vertical','Matrix'},'Vertical or matrix form?',1);
            else
                CompileFormat='Matrix';
            end
            ndim=1;
        else
            ndim=2;
            CompileFormat='Matrix';
        end
    end        
end

end

