function [SLinds,numVoxels,OutCoords] = Mask2Searchlight(Mask,SLrad,SLThresh,Shape,ThreshMask)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 2
        Shape='Sphere';
        SLThresh=5;
        ThreshMask=Mask;
    end
    if nargin == 3
        Shape='Sphere';
        ThreshMask=Mask;
    end
    if nargin == 4
        ThreshMask=Mask;
    end    
    coords=single(mat2coords(Mask));
    numSLs=size(coords,1);
    %generate searchlight shape (spherical)
    SLdiam=(SLrad*2)+1;
    if strcmpi(Shape,'Sphere')
        SLmat=zeros(SLdiam,SLdiam,SLdiam);
        centercoord=[SLrad+1,SLrad+1,SLrad+1]';
        for inix=1:SLdiam
            for iniy=1:SLdiam
                for iniz = 1:SLdiam
                    tstcoord=[inix,iniy,iniz];
                    if dist(tstcoord,centercoord)<=SLrad
                        SLmat(inix,iniy,iniz)=1;
                    end
                end
            end
        end
    else
        SLmat=ones(SLdiam,SLdiam,SLdiam);
    end
    SLcoords=single(mat2coords(SLmat))-(SLrad+1);
    SLsize=size(SLcoords,1);
    %generate searchlight matrix
    SLinds=cell(numSLs,1);
    numVoxels=zeros(numSLs,1,'single');
    removeInd=numVoxels;
    OutCoords=zeros(numSLs,3,'single');
    for i = 1:numSLs
        tempCoords=SLcoords+repmat(coords(i,:),[SLsize,1]);
        tempCoords(min(tempCoords,[],2)<1,:)=[];
        maxFilter=sum(single((repmat(size(Mask),[size(tempCoords,1),1])-tempCoords)<0),2)~=0;
        tempCoords(maxFilter,:)=[];
        try
            tempMask=coords2mat(tempCoords,Mask,ones(size(tempCoords,1),1));        
        catch
            removeInd(i,1)=1;
            disp(i)
            continue
        end
        
        tempInds=single(find(tempMask.*ThreshMask));
        numCoords=length(tempInds);
        if numCoords > SLThresh            
            SLinds{i,1}=tempInds;
            numVoxels(i,1)=length(SLinds{i,1});
            OutCoords(i,:)=coords(i,:);
        else
            removeInd(i,1)=1;
        end     
    end
    SLinds(removeInd==1,:)=[];
    numVoxels(removeInd==1,:)=[];
    OutCoords(removeInd==1,:)=[];
end
