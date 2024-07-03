function vectormat = mat2uppertriuvectormat(mat)
%Written by David Rothlein
if istable(mat)
    labelPairs1 = labels2uppertriuvectorlabels( mat.Properties.VariableNames');
    mat=table2array(mat);
    indMat=ones(size(mat,1),size(mat,1));
    indMat=triu(indMat,1);
    vectormat=zeros(size(indMat(indMat~=0),1),size(mat,3));
    for i = 1:size(mat,3)
        tempmat=triu(mat(:,:,i),1);
        vectormat(:,i)=tempmat(indMat~=0);
    end
    vectormat=array2table(vectormat','VariableNames',labelPairs1);
else
    indMat=ones(size(mat,1),size(mat,1));
    indMat=triu(indMat,1);
    vectormat=zeros(size(indMat(indMat~=0),1),size(mat,3));
    
    for i = 1:size(mat,3)
        tempmat=triu(mat(:,:,i),1);
        vectormat(:,i)=tempmat(indMat~=0);
    end
end    