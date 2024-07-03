function [y] = tertile(x)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p=quantile(x,2);
y=x;
y(x<=p(1,1))=1;
y(x>p(1,1)&x<=p(1,2))=2;
y(x>p(1,2))=3;
end

