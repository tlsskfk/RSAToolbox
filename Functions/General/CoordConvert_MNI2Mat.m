function [MatCoords] = CoordConvert_MNI2Mat(MNICoords,MatSize)
%Written by David Rothlein
%Linear translation from MNI coordinates to Matlab matrix indicies

if nargin==1
    MatSize=[61,73,61];
end
% MAT [1,1,1] = MNI [-91,-126,-72]
% MAT [61,73,61] = MNI [90,91,109]
MatCoords=repmat(MatSize,[size(MNICoords,1),1]);
MNIBox=[-91,-126,-72;90,91,109]; 
MNIsize=range(MNIBox);
MNICoords = MNICoords + repmat(abs(MNIBox(1,:)),[size(MNICoords,1),1]); %translate to positive coordinates
MNICoords = MNICoords./repmat(MNIsize,[size(MNICoords,1),1]); % rescale to be fraction of total range
MatCoords=round(MatCoords.*MNICoords);
MatCoords(MatCoords(:,1)>MatSize(1,1),1)=MatSize(1,1);
MatCoords(MatCoords(:,2)>MatSize(1,2),2)=MatSize(1,2);
MatCoords(MatCoords(:,3)>MatSize(1,3),3)=MatSize(1,3);
MatCoords(MatCoords<1)=1;
   
end

