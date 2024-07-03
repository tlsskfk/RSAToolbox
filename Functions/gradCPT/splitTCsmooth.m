function [ smoothTC ] = splitTCsmooth( TC,replaceNum )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    replaceNum=0;
end

if TC(1,1)==replaceNum
    TC(1,1)=randi(2,1);
end

smoothTC=TC;
ind=find(TC==replaceNum);

for i = ind'
   smoothTC(i,1)=smoothTC(i-1,1);
end

