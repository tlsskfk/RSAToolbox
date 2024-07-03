function [Prefix,Suffix] = StringMatch(StringCell,splitChar)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

StringCell2=StringCell;
Inds=strfind(StringCell,splitChar);
Prefix=[];
Suffix=[];
i=1;
IndSizes=cellfun(@length,Inds);
x=find(IndSizes==max(IndSizes));

while i == 1
    if ~ isempty(Inds{i,1})
        for j = 1:length(Inds{i,1})                
            if j == 1
                testStr=StringCell{i,1}(1,1:Inds{i,1}(1,j));
                strReps=contains(StringCell,testStr);
                numStrReps=sum(strReps,1);
                if length(Inds{i,1})==1
                    Prefix =[Prefix;{testStr}];
                    StringCell(strReps,:)=[];
                    break   
                end
            elseif numStrReps == sum(contains(StringCell,StringCell{i,1}(1,1:Inds{i,1}(1,j))),1)
                testStr=StringCell{i,1}(1,1:Inds{i,1}(1,j));
                if j == length(Inds{i,1})
                    Prefix =[Prefix;{testStr}];
                    StringCell(strReps,:)=[];
                end
            else
                Prefix =[Prefix;{testStr}];
                StringCell(strReps,:)=[];
                break
            end

        end 
    else
        testStr=StringCell{i,1};
        strReps=contains(StringCell,testStr);
        Prefix =[Prefix;{testStr}];
        StringCell(strReps,:)=[];        
    end
    if isempty(StringCell)
        i=2;
    else
        Inds=strfind(StringCell,splitChar);
    end
end

for i = 1:size(Prefix,1)
    StringCell2=strrepCell(StringCell2,Prefix{i,1},''); 
end

Suffix=unique(StringCell2);

end
