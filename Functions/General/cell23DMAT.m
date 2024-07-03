function [ mat ] = cell23DMAT( cell1 )
%Written by David Rothlein
a=size(cell1,1);
b=size(cell1,2);
matsize=size(cell1{1,1});
if a>b
    mat=zeros(matsize(1,1),matsize(1,2),a);
    for i = 1:a
        mat(:,:,i)=cell1{i,1};
    end;
else
    mat=zeros(matsize(1,1),matsize(1,2),b);
    for i = 1:b
        mat(:,:,i)=cell1{1,i};
    end;    
end

