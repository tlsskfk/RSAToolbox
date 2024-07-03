function [ bins ] = makebins( n,numBins )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bins =cell(3,numBins);
bin=(n)/numBins;
allbins=round(1:bin:n);

for i = 1:numBins
    vec = zeros(n,1);
    if i == numBins
        bins{1,i}=[allbins(1,i),n];
        vec(bins{1,i}(1,1):bins{1,i}(1,2),:)=1;
        bins{2,i}=vec;
        bins{3,i}=sum(vec(:));
    else    
        bins{1,i}=[allbins(1,i),allbins(1,i+1)-1];
        vec(bins{1,i}(1,1):bins{1,i}(1,2),:)=1;
        bins{2,i}=vec;
        bins{3,i}=sum(vec(:));
    end
end
