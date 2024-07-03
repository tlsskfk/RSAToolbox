function [ spherecoords ] = coord2spherecoords( coord,rad )
%Written by David Rothlein
diam=(rad*2)+1;
mat=zeros(diam,diam,diam);
centercoord=[rad+1,rad+1,rad+1]';
for inix=1:diam
    for iniy=1:diam
        for iniz = 1:diam
            tstcoord=[inix,iniy,iniz];
            if dist(tstcoord,centercoord)<=rad
                mat(inix,iniy,iniz,:)=1;
            end;
        end;
    end;
end;

%coordmat=zeros((coord(1,1)+rad),(coord(1,2)+rad),(coord(1,3)+rad);
[sphere,~]=mat2coords(mat);
spherecoords=sphere;
spherecoords(:,1)=sphere(:,1)+(coord(1,1)-(1+rad));
spherecoords(:,2)=sphere(:,2)+(coord(1,2)-(1+rad));
spherecoords(:,3)=sphere(:,3)+(coord(1,3)-(1+rad));




