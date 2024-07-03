function [SimMat] = GroupVector2SimMat(GroupVector)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
GroupVector=GroupVector(:);
SimMat=single(repmat(GroupVector,[1,length(GroupVector)])==repmat(GroupVector',[length(GroupVector),1]));

end

