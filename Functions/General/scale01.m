function [ scaled_vector ] = scale01( vector )
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
scaled_vector = (vector - min(vector)) / ( max(vector) - min(vector) );

end

