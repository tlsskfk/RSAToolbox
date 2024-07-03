function [ skip ] = AppendSave( SaveName,VarName,SaveVar,SaveType,AppendDim )

%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    AppendDim = ndims(SaveVar)+1;
    SaveType='Overwrite';
end
if nargin == 4
    AppendDim = ndims(SaveVar)+1;
end

skip=0;
if exist(SaveName,'file')
    if strcmpi(SaveType,'Append')
        try
            temp=load(SaveName,VarName);
            SaveVar=cat(AppendDim,SaveVar,temp.(VarName));
        catch
            disp(['Append_Fail: ',SaveName,' ',VarName]);
            SaveName=[SaveName,'_AF'];
        end
    elseif strcmpi(SaveType,'Skip')
        skip=1;
    end
end
        
if skip==0
    eval([VarName,'=SaveVar;']);
    try
        save(SaveName,VarName,'-append');
    catch
        save(SaveName,VarName,'-append','-v7.3')
    end
end


end

