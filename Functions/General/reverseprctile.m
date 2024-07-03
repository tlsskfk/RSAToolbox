function [perctiles] = reverseprctile(vector,vals)
%Written by David Rothlein
vals=vals(:)';
if any(size(vector)==1)
    if size(vector,1)==1
        vector=vector(:);
    end
    ShareVectors = 1;
else
    if any(size(vector)==length(vals))       
        if size(vector,2)~=length(vals) 
            vector=vector';
        end 
        ShareVectors=0;
    else
        vector=vector(:);
        ShareVectors = 1;
    end
end

perctiles=zeros(1,size(vals,2));
vectorsize=size(vector,1);
for i = 1:size(vals,2)
    val=vals(1,i);
    if ShareVectors==0
        tempvector=vector(:,i);
    else
        tempvector=vector;
    end
    tempvector=[tempvector;val];
    tempvector=sort(tempvector);
    valmatch=tempvector==val;
    perctiles(1,i)=sum((find(valmatch)-1))/(sum(valmatch,1)*vectorsize);
    if perctiles(1,i)==0
       perctiles(1,i)=1/vectorsize;
    elseif perctiles(1,i)==1
       perctiles(1,i)=(vectorsize-1)/vectorsize;
    end
end
%zvals=icdf('normal',perctiles,0,1);