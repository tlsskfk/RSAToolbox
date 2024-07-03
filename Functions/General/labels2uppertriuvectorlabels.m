function [ labelPairs1,labelPairs2,labelPairsCell,nodes1,nodes2 ] = labels2uppertriuvectorlabels( labels )
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numLabels=length(labels);
tempmat=ones(length(labels));
tempmat=triu(tempmat,1);

labelPairs1=cell((numLabels^2-numLabels)/2,1);
labelPairs2=cell((numLabels^2-numLabels)/2,1);
labelPairsCell=cell((numLabels^2-numLabels)/2,2);
if size(labels,2)==length(labels)
    labels=labels';
end

labels1=repmat(labels,[1,numLabels]);
labels1=labels1(tempmat~=0);
labels2=repmat(labels',[numLabels,1]);
labels2=labels2(tempmat~=0);

for i = 1:length(labels1)
    labelPairs1{i,1}=[labels2{i,1},'_2_',labels1{i,1}];
%    labelPairs1{i,1}=strrep(strjoin(labelPairs1{i,1}),' ','');
    labelPairs2{i,1}=[labels1{i,1},'_2_',labels2{i,1}];
 %   labelPairs2{i,1}=strrep(strjoin(labelPairs2{i,1}),' ','');
    labelPairsCell(i,:)={labels2{i,1},labels1{i,1}};
end
nodes1=nan(size(labelPairsCell,1),2);
nodes2= nodes1;
if length(labels)<=length(labelPairs2)
    for i = 1:length(labels)
       nodes1(ismember(labels2,labels{i,1}),1)=i;
       nodes1(ismember(labels1,labels{i,1}),2)=i;
       nodes2(ismember(labels1,labels{i,1}),1)=i;
       nodes2(ismember(labels2,labels{i,1}),2)=i;   
    end    
else
    nodes1=[2,1];
    nodes2=[1,2];
end

end
