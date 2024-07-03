function mat = coords2mat(coords,matrix,vals)
%Written by David Rothlein
matrix=matrix*0;
numMats=size(vals,2);
numDims=size(coords,2);

if numDims==2
    i=sub2ind(size(matrix),coords(:,1),coords(:,2));
elseif numDims==3
    i=sub2ind(size(matrix),coords(:,1),coords(:,2),coords(:,3));
elseif numDims==4
    i=sub2ind(size(matrix),coords(:,1),coords(:,2),coords(:,3),coords(:,4));
elseif numDims==5
    i=sub2ind(size(matrix),coords(:,1),coords(:,2),coords(:,3),coords(:,4),coords(:,5));
elseif numDims==6
    i=sub2ind(size(matrix),coords(:,1),coords(:,2),coords(:,3),coords(:,4),coords(:,5),coords(:,6));
elseif numDims==7
    i=sub2ind(size(matrix),coords(:,1),coords(:,2),coords(:,3),coords(:,4),coords(:,5),coords(:,6),coords(:,7));
end

if numMats>1
    mat = cell(1,numMats);
    for n = 1:numMats
        TempMat=matrix;
        TempMat(i)=vals(:,n);
        mat{1,n}=TempMat;
    end
    mat = cell2nDMAT(mat);
else
    mat=matrix;
    mat(i)=vals;
end


    