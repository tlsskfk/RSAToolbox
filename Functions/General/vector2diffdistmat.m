function [diffpropmat] = vector2diffdistmat( v )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
diffmat=abs(bsxfun(@minus,v,v'));
maxv=max(diffmat(:));
diffpropmat=1-(diffmat/maxv);
end

