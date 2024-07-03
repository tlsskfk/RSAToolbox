function [ ellipsoidcoords ] = coord2ellipsoidcoords( coord,radx,rady,radz )
%Written by David Rothlein
diamx=(radx*2)+1;
diamy=(rady*2)+1;
diamz=(radz*2)+1;
mat=zeros(diamx,diamy,diamz);
centercoord=[radx+1,rady+1,radz+1];
for inix=1:diamx
    for iniy=1:diamy
        for iniz = 1:diamz
            if abs(inix-centercoord(1,1))<=radx || abs(iniy-centercoord(1,2))<=rady || abs(iniz-centercoord(1,3))<=radz
                mat(inix,iniy,iniz,:)=1;
            end;
        end;
    end;
end;

%coordmat=zeros((coord(1,1)+rad),(coord(1,2)+rad),(coord(1,3)+rad);
[sphere,~]=mat2coords(mat);
ellipsoidcoords=sphere;
ellipsoidcoords(:,1)=sphere(:,1)+(coord(1,1)-(1+radx));
ellipsoidcoords(:,2)=sphere(:,2)+(coord(1,2)-(1+rady));
ellipsoidcoords(:,3)=sphere(:,3)+(coord(1,3)-(1+radz));




