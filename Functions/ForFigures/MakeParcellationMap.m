function [ useMap ] = MakeParcellationMap( Vals,ParcellationMask)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numMaps=size(Vals,2);
numNodes=size(Vals,1);

useMap=repmat(ParcellationMask*0,[1,1,1,numMaps]);

for mapNum = 1:numMaps
    tempMap=ParcellationMask*0;
    for nodeNum = 1:numNodes
        tempMap(ParcellationMask==nodeNum)=Vals(nodeNum,mapNum);
    end
    useMap(:,:,:,mapNum)=tempMap;
end

end

