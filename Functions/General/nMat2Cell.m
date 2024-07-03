function MatCell=nMat2Cell(mat,nDim)
%Written by David Rothlein
if nargin == 1
    nDim='';
end;
if isnumeric(nDim)==0 || isempty(nDim)==1
    nDim=length(size(mat));
end;
totalDims=length(size(mat));
MatCell=cell(size(mat,nDim),1);

for i = 1:size(mat,nDim)
    if totalDims==2
        if nDim==1
            MatCell{i,1}=mat(i,:);
        elseif nDim==2
            MatCell{i,1}=mat(:,i);
        end;
    elseif totalDims==3
        if nDim==1
            MatCell{i,1}=mat(i,:,:);
        elseif nDim==2
            MatCell{i,1}=mat(:,i,:);
        elseif nDim==3
            MatCell{i,1}=mat(:,:,i);
        end;
    elseif totalDims==4
        if nDim==1
            MatCell{i,1}=mat(i,:,:,:);
        elseif nDim==2
            MatCell{i,1}=mat(:,i,:,:);
        elseif nDim==3
            MatCell{i,1}=mat(:,:,i,:);
        elseif nDim==4
            MatCell{i,1}=mat(:,:,:,i);            
        end;
    elseif totalDims==5
        if nDim==1
            MatCell{i,1}=mat(i,:,:,:,:);
        elseif nDim==2
            MatCell{i,1}=mat(:,i,:,:,:);
        elseif nDim==3
            MatCell{i,1}=mat(:,:,i,:,:);
        elseif nDim==4
            MatCell{i,1}=mat(:,:,:,i,:);  
        elseif nDim==5
            MatCell{i,1}=mat(:,:,:,:,i); 
        end;
    elseif totalDims==6
        if nDim==1
            MatCell{i,1}=mat(i,:,:,:,:,:);
        elseif nDim==2
            MatCell{i,1}=mat(:,i,:,:,:,:);
        elseif nDim==3
            MatCell{i,1}=mat(:,:,i,:,:,:);
        elseif nDim==4
            MatCell{i,1}=mat(:,:,:,i,:,:);  
        elseif nDim==5
            MatCell{i,1}=mat(:,:,:,:,i,:); 
        elseif nDim==6
            MatCell{i,1}=mat(:,:,:,:,:,i); 
        end;
    end;
end
