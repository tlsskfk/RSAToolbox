function [ r,p,n ] = nancorr( x,y,type )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    type='Pearson';
end;
numx=size(x,2);
numy=size(y,2);

r=zeros(numx,numy);
p=zeros(numx,numy);
n = zeros(numx,numy);
count=1;
for i = 1:numx
    a=x(:,i);
    temp_ind1=isnan(a);
    
    for j = 1:numy
        b = y(:,j);
        temp_ind2 = isnan(b);
        temp_ind = temp_ind1+temp_ind2;
        temp_ind = temp_ind == 0;
        a2=a(temp_ind);
        b2 = b(temp_ind);
        if isempty(a2) || length(a2)==1
           r(i,j)=nan;
           p(i,j)=nan;
        else
            [r(i,j),p(i,j)]=corr(a2,b2,'type',type);
            n(i,j) = size(a2,1);
        %count/total
        end
        count=count+1;
        
    end;
        

end

