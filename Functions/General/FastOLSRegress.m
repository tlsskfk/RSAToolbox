function [T,B,SE,E,DFE] = FastOLSRegress(Y,X)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
idX=eye(size(X,2));
idY=eye(size(Y,2));
DFE=size(X,1)-size(X,2);
B=Y'/X';
E=(Y-X*B')'*(Y-X*B')/DFE;
E=E(idY==1);
SE=idX*inv(X'*X)*idX;
SE=sqrt(SE(idX==1)*E')';
T=B./SE;
end

