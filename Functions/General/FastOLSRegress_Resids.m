function [Resids] = FastOLSRegress_Resids(Y,X)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numY=size(Y,2);
Resids=single(Y*0);
parfor i = 1:numY 
    tempY=Y(:,i)
    B=tempY'/X';
    Resids(:,i)=tempY-[B*X']';
    warning('off','all'); 
end
end
