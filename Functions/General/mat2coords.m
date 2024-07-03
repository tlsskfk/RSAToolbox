function [coords,vals] = mat2coords(matrix)
%Written by David Rothlein
matrix(isnan(matrix)==1)=0;
matrix(isinf(matrix)==1)=0;
i=find(matrix);
vals=matrix(i);
vals=single(vals);
if length(size(matrix))==2
    [x,y]=ind2sub(size(matrix),i);
    coords=[x,y];
elseif length(size(matrix))==3
    [x,y,z]=ind2sub(size(matrix),i);
    coords=[x,y,z];
elseif length(size(matrix))==4
    [x,y,z,a]=ind2sub(size(matrix),i);
    coords=[x,y,z,a];
elseif length(size(matrix))==5
    [x,y,z,a,b]=ind2sub(size(matrix),i);
    coords=[x,y,z,a,b];
elseif length(size(matrix))==6
    [x,y,z,a,b,c]=ind2sub(size(matrix),i);
    coords=[x,y,z,a,b,c];
elseif length(size(matrix))==7
    [x,y,z,a,b,c,d]=ind2sub(size(matrix),i);
    coords=[x,y,z,a,b,c,d];
end

if max(coords(:))<256
    coords=uint8(coords);
elseif max(coords(:))<65535
    coords=uint16(coords);
else
    coords=uint32(coords);
end

