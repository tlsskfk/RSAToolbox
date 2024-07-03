function [ Mat ] = CellConvert( CellData )
%CELLCONVERT Summary of this function goes here
%   Detailed explanation goes here
Mat=nan(length(CellData));

for i = 1:size(CellData,1)
    for j= i+1:size(CellData,2)
        if isnumeric(CellData{i,j})
            Mat(i,j)=CellData{i,j};
            Mat(j,i)=CellData{i,j};
        end
    end
end

