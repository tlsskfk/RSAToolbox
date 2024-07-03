function [newcoords,ifilter]=filtercoords(coords,max_x,max_y,max_z,min)
%Written by David Rothlein
filterx=coords(:,1)>max_x;
filtery=coords(:,2)>max_y;
filterz=coords(:,3)>max_z;
filtermx=coords(:,1)<min;
filtermy=coords(:,2)<min;
filtermz=coords(:,3)<min;
filter=filterx+filtery+filterz+filtermx+filtermy+filtermz;
ifilter=filter>0;
coords(filter>0,:)=[];
newcoords=coords;