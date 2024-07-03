function NewMap = SaveBrik_3mmMNI(mat,labels,SaveName)
%Written by David Rothlein
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if iscell(mat)
     [ mat ] = cell2nDMAT( mat );
end

if isempty(labels)
    labels=cell(size(mat,4),1);
    for i = 1:size(mat,4)
        labels{i,1}=num2str(i);
    end
end

mat=single(mat);
if ~iscell(labels)
    labels={labels};
end
labels=labels(:);
numLabels=length(labels);
BrikLabels=[];

for i =1:numLabels
    BrikLabels=[BrikLabels,labels{i,1}];
    if i<numLabels
        BrikLabels=[BrikLabels,'~'];
    end
end

LoadOpt.Format = 'matrix';
BrikName='Parcellations/RawFiles/master_3mmMNI';
[~, ~, Info, ~] = BrikLoad (BrikName, LoadOpt);
Info.BRICK_LABS=BrikLabels;
SaveOpt.Prefix=SaveName;
SaveOpt.OverWrite='y';
SaveOpt.View='+tlrc';


WriteBrik (mat, Info, SaveOpt);
