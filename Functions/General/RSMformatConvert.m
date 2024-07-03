function [ RSMvert,RSMmat ] = RSMformatConvert( Data )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if iscell(Data)
    Data=squeeze(cell2nDMAT(Data));
end
if ndims(Data)==4
    numROIs=size(Data,4);
    numSS=size(Data,3);
    RSMsize=size(Data,2);
    RSMvert=zeros((RSMsize^2-RSMsize)/2,numSS,numROIs);
    RSMs=Data;
    for i = 1:numROIs
        for j = 1:numSS
            RSMvert(:,j,i)=mat2uppertriuvectormat(Data(:,:,j,i));
            RSMs(:,:,j,i)=triRSMs2SymRSMs(RSMs(:,:,j,i));
        end
    end
elseif ndims(Data)==3 && (size(Data,1)~=size(Data,2) || (~isnan(Data(1,1)) && Data(1,1)~=1))
    RSMvert=Data;
    numROIs=size(Data,3);
    vRSMsize=size(Data,1);
    numSS=size(Data,2);
    RSMsize=roots([1,-1,-(2*vRSMsize)]);
    RSMsize(RSMsize<0,:)=[];
    RSMs=zeros(RSMsize,RSMsize,numSS,numROIs);
    for i = 1:numROIs
        for j = 1:numSS
            [ RSMs(:,:,j,i) ] = vertRSM2SymRSM( RSMvert(:,j,i) );
        end    
    end
elseif ndims(Data)==3
    numSS=size(Data,3);
    numROIs=1;
    RSMsize=size(Data,2);
    RSMvert=zeros((RSMsize^2-RSMsize)/2,numSS);
    RSMs=Data;
    for j = 1:numSS
        RSMvert(:,j)=mat2uppertriuvectormat(Data(:,:,j));
        RSMs(:,:,j)=triRSMs2SymRSMs(RSMs(:,:,j));
    end  
elseif ndims(Data)>=2
    RSMvert=Data;
    numROIs=1;
    vRSMsize=size(Data,1);
    numSS=size(Data,2);
    RSMsize=roots([1,-1,-(2*vRSMsize)]);
    RSMsize(RSMsize<0,:)=[];
    RSMs=zeros(RSMsize,RSMsize,numSS);
    for j = 1:numSS
        [ RSMs(:,:,j) ] = vertRSM2SymRSM( RSMvert(:,j) );
    end       
end
RSMmat=RSMs;
end

