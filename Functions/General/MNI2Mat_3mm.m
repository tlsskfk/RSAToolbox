function [MatCoords] = MNI2Mat(MNICoords,MatSize)
%Written by David Rothlein
% MAT [1,1,1] = MNI [-91,-126,-72]
% MAT [61,73,61] = MNI [90,91,109]

    MatCoords=repmat([61,73,61],[size(MNICoords,1),1]);
    MNIBox=[-91,-126,-72;90,91,109]; 
    MNIsize=range(MNIBox);
    MNICoords = MNICoords + repmat(abs(MNIBox(1,:)),[size(MNICoords,1),1]); %translate to positive coordinates
    MNICoords = MNICoords./repmat(MNIsize,[size(MNICoords,1),1]); % rescale to be fraction of total range
    %MNICoords(:,2)=1-MNICoords(:,2); %Unflip Y dimension
    MatCoords=round(MatCoords.*MNICoords) +1;
    MatCoords(MatCoords(:,1)>61,1)=61;
    MatCoords(MatCoords(:,2)>73,2)=73;
    MatCoords(MatCoords(:,3)>61,3)=61;
    
end

