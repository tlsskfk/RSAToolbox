function [ X ] = scaleVals( X,dim,lo,hi)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    dim=1;
    lo=0;
    hi=1;
elseif nargin == 2
    lo=0;
    hi=1;
elseif nargin==3
    hi=1;
end

if dim==1
    minX=repmat(min(X),[size(X,1),1]);
    maxX=repmat(max(X),[size(X,1),1]);
    rangeX=maxX-minX;
    X=(((X-minX)./rangeX)*(hi-lo))+lo;    
elseif dim == 2
    X=X';
    minX=repmat(min(X),[size(X,1),1]);
    maxX=repmat(max(X),[size(X,1),1]);
    rangeX=maxX-minX;
    X=(((X-minX)./rangeX)*(hi-lo))+lo; 
    X=X';
else
    minX=min(X(:));
    maxX=max(X(:));
    rangeX=maxX-minX;
    X=(((X-minX)/rangeX)*(hi-lo))+lo;
end
    
%X(isnan(X))=lo+((hi-lo)/2);
end

