function [backtime] = ttt2backtime(ttt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x=size(ttt,1)/4;
y=size(ttt,2)*4;
backtime=[];
tempVec=[];
count=1;
for i = 1:size(ttt,1)
    if count==1
        tempVec=ttt(i,:);
        count=count+1;
    elseif count==4
        tempVec=[tempVec,ttt(i,:)];
        backtime=[backtime;tempVec];
        count=1;
    else
        tempVec=[tempVec,ttt(i,:)];
        count=count+1;
    end
end    

end