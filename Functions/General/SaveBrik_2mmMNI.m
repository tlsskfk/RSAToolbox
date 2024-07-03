function NewMap = SaveBrik_2mmMNI(mat,labels,SaveName)
%Written by David Rothlein
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if iscell(mat)
     [ mat ] = cell2nDMAT( mat );
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
BrikName='Parcellations/RawFiles/master_2mmMNI+tlrc';
[~, ~, Info, ~] = BrikLoad (BrikName, LoadOpt);
Info.BRICK_LABS=BrikLabels;
SaveOpt.Prefix=SaveName;
SaveOpt.OverWrite='y';
SaveOpt.View='+tlrc';


WriteBrik (mat, Info, SaveOpt);
