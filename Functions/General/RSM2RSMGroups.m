function [wRSMs,xRSMs] = RSM2RSMGroups(RSMs,vertIn)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin ==1
    vertIn=1;
end

if vertIn == 1
    numCond=[];
    numROIs=size(RSMs,2);
    numSS=size(RSMs,3);
    wRSMs=[];
    xRSMs=[];
    for i = 1:numSS
        [ SymRSMs ] = vertRSM2SymRSM( RSMs(:,:,i),numCond,1 );
        [ wRSM,xRSM,Diags] = RSMParcer( SymRSMs,2 );
        wRSM=reshape(wRSM,[size(wRSM,1)*size(wRSM,2),size(wRSM,3)]);
        xRSM=[reshape(xRSM,[size(xRSM,1)*size(xRSM,2),size(xRSM,3)]);reshape(Diags,[size(Diags,1)*size(Diags,2),size(Diags,3)])];
        if isempty(wRSMs)
            numCond=size(SymRSMs,1);
            wRSMs=zeros(size(wRSM,1),numROIs,numSS,'single');
            xRSMs=zeros(size(xRSM,1),numROIs,numSS,'single');
        end
        wRSMs(:,:,i)=wRSM;
        xRSMs(:,:,i)=xRSM;
    end
else
    numROIs=size(RSMs,3);
    numSS=size(RSMs,4);
    wRSMs=[];
    xRSMs=[];
    numCond=size(RSMs,1);
    for i = 1:numSS
        [ wRSM,xRSM,Diags] = RSMParcer( RSMs(:,:,:,i),2 );
        wRSM=reshape(wRSM,[size(wRSM,1)*size(wRSM,2),size(wRSM,3)]);
        xRSM=[reshape(xRSM,[size(xRSM,1)*size(xRSM,2),size(xRSM,3)]);reshape(Diags,[size(Diags,1)*size(Diags,2),size(Diags,3)])];
        if isempty(wRSMs)
            wRSMs=zeros(size(wRSM,1),numROIs,numSS,'single');
            xRSMs=zeros(size(xRSM,1),numROIs,numSS,'single');
        end
        wRSMs(:,:,i)=wRSM;
        xRSMs(:,:,i)=xRSM;
    end
end

end

