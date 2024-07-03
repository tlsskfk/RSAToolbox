function [IndvPs,PairwisePs] = AddSigMarkers(IndvPs,PairwisePs,XsForStats,MaxPts)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

YUnit=range(ylim)/100;
if ~isempty(IndvPs)
    [ ~,IndvPs ] = SigSeg(IndvPs);
    %Plot Indv astricies
    for i = 1:size(XsForStats,2)
        for j = 1:size(XsForStats,1)
            if ~isempty(IndvPs{j,i})
                UseX=XsForStats(j,i);
                UseY=MaxPts(j,i)+YUnit*2;
                text(double(UseX),double(UseY),IndvPs{j,i},'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Interpreter','none');
                hold on    
            end
        end
    end
end
if ~isempty(PairwisePs)        
    if ~isempty(PairwisePs)
        for i = 1:size(XsForStats,2)
            [~,tempPairwisePs] = SigSeg(PairwisePs(:,:,i));
            for j = 1:size(tempPairwisePs,1)
                UseX1=XsForStats(j,i);
                if j ~= size(tempPairwisePs,1)
                    for k = j+1:size(XsForStats,1)
                        if ~isempty(tempPairwisePs{j,k})
                            UseX2=XsForStats(k,i);                        
                            UseX=(UseX1+UseX2)/2;
                            UseYLine=max(MaxPts(:,i))+YUnit*(3+((k-j)*3)+j);
                            UseYSig = UseYLine+YUnit;
                            text(double(UseX),double(UseYSig),tempPairwisePs{j,k},'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Interpreter','none');
                            hold on    
                            plot([UseX1,UseX2],[UseYLine,UseYLine],'k');
                            hold on
                        end
                    end
                end
            end
        end    
    end    
end
end

function [ SigMat,SigLabelMat ] = SigSeg( pMat,pThreshs,signType,ThreshLabels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin==1
        pThreshs=[0.05,0.01,0.005,0.001];
        signType='pos';
        for i = 1:length(pThreshs)
            ThreshLabels{1,i}=repmat('*',[1,i]);
        end
    end
    if nargin==2
        signType='pos';
        for i = 1:length(pThreshs)
            ThreshLabels{1,i}=repmat('*',[1,i]);
        end
    end
    if nargin==3
        for i = 1:length(pThreshs)
            ThreshLabels{1,i}=repmat('*',[1,i]);
        end
    end

    SigMat=pMat*0;
    for i = 1:length(pThreshs)
        if strcmpi(signType,'both')
            SigMat=SigMat+single(pMat<=pThreshs(1,i)|pMat>=1-pThreshs(1,i));
        elseif strcmpi(signType,'pos')
            SigMat=SigMat+single(pMat<=pThreshs(1,i));
        elseif strcmpi(signType,'neg')
            SigMat=SigMat+single(pMat>=1-pThreshs(1,i));
        end
    end

    SigLabelMat=cell(size(SigMat));
    for i =1:size(SigMat,1)
        for j = 1:size(SigMat,2)
            if SigMat(i,j)>0
                SigLabelMat{i,j}=ThreshLabels{1,SigMat(i,j)};
            else
                SigLabelMat{i,j}=[];
            end
        end
    end
end

