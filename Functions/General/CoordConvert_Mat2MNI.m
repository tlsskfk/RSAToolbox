function [MNICoords] = CoordConvert_Mat2MNI(MatCoords,MatSize)
%Written by David Rothlein
%Linear translation from MNI coordinates to Matlab matrix indicies
if nargin==1
    MatSize=[61,73,61];
end
% MAT [1,1,1] = MNI [-91,-126,-72]
% MAT [61,73,61] = MNI [90,91,109]
% MatCoords=MatCoords-1;
% MatSize=MatSize-1;
MNIBox=[-91,-126,-72;90,91,109]; 
MNIsize=range(MNIBox);
MNICoords=repmat(MNIsize,[size(MatCoords,1),1]);
MatCoords = MatCoords./repmat(MatSize,[size(MatCoords,1),1]); %  rescale matcoords to be fraction of total range
MNICoords=round(MatCoords.*MNICoords);
MNICoords = MNICoords + repmat(MNIBox(1,:),[size(MNICoords,1),1]); %shift to MNI coordinates

MNICoords(MNICoords(:,1)<MNIBox(1,1),1)=MNIBox(1,1);
MNICoords(MNICoords(:,2)<MNIBox(1,2),2)=MNIBox(1,2);
MNICoords(MNICoords(:,3)<MNIBox(1,3),3)=MNIBox(1,3);
MNICoords(MNICoords(:,1)>MNIBox(2,1),1)=MNIBox(2,1);
MNICoords(MNICoords(:,2)>MNIBox(2,2),2)=MNIBox(2,2);
MNICoords(MNICoords(:,3)>MNIBox(2,3),3)=MNIBox(2,3);
    
end

