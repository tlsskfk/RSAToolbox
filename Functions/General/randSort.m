function [ SortInd ] = randSort( Vals )
%Written by David Rothlein
%Sorts but randomly assigns the rank of equal values

uVals=unique(Vals);
[Vals,SortInd]=sort(Vals);
for i = uVals(:)'
    tempInd=SortInd(Vals==i);
    tempInd=tempInd(randperm(length(tempInd)));
    SortInd(Vals==i)=tempInd;
end

