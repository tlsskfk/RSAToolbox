function [sctPlot,Xs] = TearDropPlot( inData,labels,varargin)
%Written by David Rothlein
numBins = VariableSetter( 'numBins',30,varargin);
yRes = VariableSetter( 'yRes',0.001,varargin);
xRes = VariableSetter( 'xRes',0.001,varargin);
MarkerSize = VariableSetter( 'MarkerSize',5,varargin);
colors = VariableSetter( 'colors',[0.5,0.5,0.5],varargin);
transparency = VariableSetter( 'transparency',0.025,varargin);
CI = VariableSetter( 'CI',[0,100],varargin);
sigColor = VariableSetter( 'sigColor',[0.3,0.3,0.3],varargin);
sigLine = VariableSetter( 'sigLine',[],varargin);
sigLineColor = VariableSetter( 'sigLineColor',[],varargin);
totalAreas = VariableSetter( 'totalAreas',ones(1,size(inData,2)),varargin);
        %  1 byN vector of values between 0 and 1 for each X indicating  
        %  percentage of DataPoints in a given plot relative to plot 
        %  with most datapoints. If all equal, defaults to 1. 
smoothness = VariableSetter( 'smoothness',0.3,varargin);
smoothType = VariableSetter( 'smoothType','loess',varargin);
xShift = VariableSetter( 'xShift',0,varargin);
        %  Scalar value (generally between 0 and 1). Shifts 
        %  reported default X positions. Differentiates overlapping plots.
meanMarker = VariableSetter( 'meanMarker','o',varargin);
meanMarkerColor = VariableSetter( 'meanMarkerColor',[0,0,0],varargin);
meanMarkerSize = VariableSetter( 'meanMarkerSize',20,varargin);
MeanMarkerFill = VariableSetter( 'MeanMarkerFill',1,varargin); 

if transparency <=0
    transparency=0.001;
end

Data=[];
widthData=[];
sigData=[];
sigLineData=[];
X_fix=0;
if size(inData,1)==1
    X_fix=1;
end
    
for i = 1:size(inData,2)
    if X_fix==1
        inData{2,i}=i;
    end
    tempData=inData{1,i};
    [n,edges]=histcounts(tempData,numBins);
    tempData=[min(tempData):yRes:max(tempData)]';
    bin = discretize(tempData,edges);
    tempWidth=[n(bin)]';
    

    
    widthRange=totalAreas(1,i)*0.9;
    
    
    Data4Sig=[];
    if isempty(smoothness)==0
        tempWidth=smooth(tempWidth,smoothness,smoothType); 
    end
    tempWidth=scaleVals(tempWidth,0,widthRange);
    for j=1:length(tempData)
        tempWidth2=[0:xRes:tempWidth(j,1)]';
        tempWidth2=(tempWidth2-mean(tempWidth2))+inData{2,i};
        tempData2=ones(length(tempWidth2),1)*tempData(j,1);
        widthData=[widthData;tempWidth2];
        Data=[Data;tempData2];
        Data4Sig=[Data4Sig;tempData2];
    end

    if length(CI)>0 
        p=prctile(inData{1,i},CI);
        if isempty(sigColor)==0
            if length(CI)>1
                tSig=single(Data4Sig<=p(1,1)) + single(Data4Sig>=p(1,end));
            else
                tSig=single(Data4Sig>=p(1,1));
            end
            sigData=[sigData;tSig];
        end        
        if isempty(sigLine)==0
            tSigLine=0;
            for z=length(p)
                [~,critValInd]=min(abs(Data4Sig-p(1,z)));
                critVal=Data4Sig(critValInd,1);
                tSigLine=tSigLine+single(Data4Sig==critVal);
            end
            sigLineData=[sigLineData;tSigLine];
        end
    end
end
sigAdd=1;
if length(CI)>0
    if isempty(sigData)==0
        colors=[colors;sigColor];
        sigData=sigData+sigAdd;
        sigAdd=0;
    end
    if isempty(sigLineData)==0
        colors=[colors;sigLineColor];
        sigData=sigLineData+sigAdd;
    end
    colors=colors(sigData,:);
end

widthData=widthData+xShift;
sctPlot=scatter(widthData,Data,MarkerSize,colors,'Marker','.');
sctPlot.MarkerEdgeAlpha = transparency;

hold on

for j=1:size(inData,2)
    Means(j,1)=nanmean(inData{1,j});
    Xs(j,1)=j;
end

Xs=Xs+xShift;
if ~isempty(meanMarker)
if MeanMarkerFill==1
    scatter(Xs,Means,meanMarkerSize,meanMarkerColor,'filled','Marker',meanMarker);
else
    scatter(Xs,Means,meanMarkerSize,meanMarkerColor,'Marker',meanMarker);
end
end

set(gca,'Xticklabel',[]);
set(gca,'XTick',[]);

if isempty(labels)==0
    yPos=min(ylim)-(range(ylim)/20);
    for L=1:size(labels,2)
        text(Xs(L,1),yPos,labels{1,L},'HorizontalAlignment','center');
    end
end