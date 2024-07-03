function [ mat ] = cell2nDMAT( cell1 )
%Written by David Rothlein
ndim=length(size(cell1{1,1}));
a=size(cell1,1);
b=size(cell1,2);
mat=[];
if a>b
    for i = 1:a
        mat=cat(ndim+1,mat,cell1{i,1});
    end
else
    for i = 1:b
        mat=cat(ndim+1,mat,cell1{1,i});
    end
end

