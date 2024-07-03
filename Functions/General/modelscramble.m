function [newmodels,rperms_out]=modelscramble(models,type,rperms,varargin)
%Written by David Rothlein
vertIn = VariableSetter( 'vertIn',0,varargin);

if vertIn == 1
    [ models ] = vertRSM2SymRSM(models);
else
    models=triRSMs2SymRSMs(models);
end
if nargin==2
    rperms=[];
end

numofmod=size(models,3);
newmodels=zeros(size(models,1),size(models,2),numofmod);
if strcmp(type,'cell')==1
    for i = 1:numofmod
        model=models(:,:,i);
        ind=model~=99999;
        modelsize=sum(ind(:));
        vals=model(ind);
        newvals=vals(randperm(modelsize)',1);
        model(model~=99999)=newvals;
        newmodels(:,:,i)=model;
    end
elseif strcmp(type,'VH')==1
    rperms_out=zeros(numofmod,size(models,1));
    for i = 1:numofmod
        model=models(:,:,i);
        if isempty(rperms)
            a=randperm(size(model,1));
        else
            a=rperms(i,:);
        end
        model=model(a,:);
        model=model(:,a);
        newmodels(:,:,i)=model;
        rperms_out(i,:)=a;
    end
elseif strcmp(type,'VH_once')==1
    if isempty(rperms)
        a=randperm(size(models,1));
    else
        a=rperms;
    end    
    rperms_out=a;
    for i = 1:numofmod
        model=models(:,:,i);        
        model=model(a,:);
        model=model(:,a);
        newmodels(:,:,i)=model;
    end
end

if vertIn == 1
    newmodels=mat2uppertriuvectormat(newmodels);
end
end