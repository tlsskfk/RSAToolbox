function [T,B,SE,E,DFE] = FastOLSRegressContrast(Y,X,cMat)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
idX=eye(size(X,2));
fMat=[idX;cMat];
idY=eye(size(Y,2));
DFE=size(X,1)-size(X,2);
B=Y'/X';
E=(Y-X*B')'*(Y-X*B')/DFE;
E=E(idY==1);
B=[B,B*cMat'];
SE=fMat*inv(X'*X)*fMat';
SE=sqrt(SE(eye(size(fMat,1))==1)*E')';
T=B./SE;
end