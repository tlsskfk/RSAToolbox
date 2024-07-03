function [UseAxes] = DataPointDropPlot_4App( PlotInfo)
%Written by David Rothlein
AllXs=[];
UsePs='Perm';
if iscell(PlotInfo.GroupName)
    GroupNames=PlotInfo.GroupName(:);
else
    GroupNames={PlotInfo.GroupName};
end
UseAxes=zeros(1,height(PlotInfo));
MaxPts=nan(height(PlotInfo),size(PlotInfo.Data,2));
for j =1:height(PlotInfo)
    inData=PlotInfo.Data(j,:);
    xRes=size(inData,2)/2000;    
    yRange = PlotInfo.yRange(j,:);
    yRes=abs(yRange(1,2)-yRange(1,1))/2000;
    AddLine = PlotInfo.AddLine(j,:);
    XLinePosition= PlotInfo.XLinePosition(j,:);
    Bootstrap = PlotInfo.Bootstrap(j,:);
    MarkerSize = PlotInfo.MarkerSize(j,:);
    MarkerType = PlotInfo.MarkerType{j,:};
    MarkerFill = PlotInfo.MarkerFill(j,:);
    colors = PlotInfo.colors(j,:);
    %transparency = PlotInfo.transparency(i,:);
    totalAreas = PlotInfo.totalAreas{j,:};
    maxArea = PlotInfo.maxArea(j,:);
    BootstrapMaxArea = PlotInfo.BootstrapMaxArea(j,:);
    bootstrap_transparency=PlotInfo.BootstrapTransparency(j,:);
    BootstrapColor = PlotInfo.BootstrapColors(j,:);
            %  1 byN vector of values between 0 and 1 for each X indicating  
            %  percentage of DataPoints in a given plot relative to plot 
            %  with most datapoints. If all equal, defaults to 1. 
    xShift = PlotInfo.xShift(j,:);
            %  Scalar value (generally between 0 and 1). Shifts 
            %  reported default X positions. Differentiates overlapping plots.
    meanMarker = PlotInfo.meanMarker{j,:};
    meanMarkerColor = PlotInfo.meanMarkerColor(j,:);
    meanMarkerSize = PlotInfo.meanMarkerSize(j,:);
    MeanMarkerFill = PlotInfo.meanMarkerFill(j,:); 
    
        
%     if transparency <=0
%         transparency=0.001;
%     end

    if strcmpi(totalAreas,'Identical')
        totalAreas=ones(1,size(inData,2));
    else
        Ns=zeros(1,size(inData,2));
        for i =1:size(inData,2)
            Ns(1,i)=length(inData{1,i});
        end
        totalAreas=max(Ns)./Ns;
    end
    BootstrapAreas=totalAreas*BootstrapMaxArea;
    totalAreas=totalAreas*maxArea;

    Data=[];
    widthData=[];
    sigData=[];
    sigLineData=[];
    X_fix=0;
    sampleDist=cell(1,size(inData,2));
    if size(inData,1)==1
        X_fix=1;
    end

    for i = 1:size(inData,2)
        if X_fix==1
            inData{2,i}=i;
        end
        tempData=inData{1,i};
        tempData(isnan(tempData),:)=[];
        tempData(isinf(tempData),:)=[];
        if isempty(tempData)
            continue
        end
        MaxPts(j,i)=max(tempData);
        numBins=round(length(tempData)/5);
        if numBins<5
            numBins=5;
        elseif numBins > 70
            numBins = 70;
        end       
        [n,edges]=histcounts(tempData,numBins);
        %tempData=[min(tempData):yRes:max(tempData)]';
        bin = discretize(tempData,edges);
        tempWidth=[n(bin)]';
        tempWidth=scale01(arrayfun(@randi,tempWidth));
        widthInd=single(tempWidth~=0);
        nonZeros=sum(widthInd(:));
        symVec=zeros(nonZeros,1)-1;
        symVec(1:round(nonZeros/2),1)=1;
        symVec(randperm(nonZeros),1)=symVec;
        widthInd(widthInd==1,1)=symVec;
        tempWidth=tempWidth.*widthInd;
        widthRange=totalAreas(1,i)*0.5;
        Data=[Data;tempData];
        widthData = [widthData;tempWidth*widthRange+inData{2,i}]; 
        if ~isempty(Bootstrap)
            bsDist=zeros(Bootstrap,1);
            for bsRep=1:Bootstrap
                bsDist(bsRep,1)=mean(tempData(randi(length(tempData),[length(tempData),1]),1));
            end
            sampleDist{1,i}=bsDist;
        end
    end

    widthData=widthData+xShift;
    if MarkerFill==1
        UseAxes(j)=scatter(widthData,Data,MarkerSize,colors,'filled','Marker',MarkerType,'DisplayName',GroupNames{j,1});
    else
        UseAxes(j)=scatter(widthData,Data,MarkerSize,colors,'Marker',MarkerType,'DisplayName',GroupNames{j,1});
    end
    %UseAxes.MarkerEdgeAlpha = transparency;
    hold on
    if ~isempty(Bootstrap)
        if max(BootstrapColor) >= 0.5
            sigBootstrapColor=BootstrapColor*0.5;
        else
            sigBootstrapColor=BootstrapColor*2;
        end
        TearDropPlot_4App(sampleDist,[],'yRes',yRes,'xRes',xRes,'transparency',bootstrap_transparency,'totalAreas',BootstrapAreas,'CI',[2.5,97.5],'meanMarker',[],'colors',BootstrapColor,'sigColor',sigBootstrapColor,'xShift',xShift);
        hold on
    end
    for k=1:size(inData,2)
        Means(k,1)=nanmean(inData{1,k});
        Xs(k,1)=inData{2,k};
    end

    Xs=Xs+xShift;
    if MeanMarkerFill==1
        scatter(Xs,Means,meanMarkerSize,meanMarkerColor,'filled','Marker',meanMarker);
        hold on
    else
        scatter(Xs,Means,meanMarkerSize,meanMarkerColor,'Marker',meanMarker);
        hold on
    end
    AllXs=[AllXs;Xs(:)'];
    if ~isempty(yRange)
        ylim(gca,yRange)
        hold on
    end
    set(gca,'Xticklabel',[]);
    hold on
    set(gca,'XTick',[]);
    hold on
    if AddLine==1
        xLine=xlim;
        plot(xLine,[XLinePosition,XLinePosition],'k');
        hold on
    end
end
LabelShare=any(PlotInfo.LabelShare==1); 
LabelStagger=max(PlotInfo.LabelStagger(:));
LabelRotate=max(PlotInfo.LabelRotate(:)) ;   
IncludeLegend=any(PlotInfo.IncludeLegend==1);  
LegendLocation=PlotInfo.LegendLocation{1,1};
XsForStats=AllXs;
if LabelShare==1
    AllXs=mean(AllXs,1);
end
Staggers=reshape(mod(1:length(AllXs(:)),LabelStagger)+1,[size(AllXs,1),size(AllXs,2)]);

for i = 1:size(AllXs,1)
    labels=PlotInfo.XLabel(i,:);
    LabelFontSize=PlotInfo.LabelFontSize(i,:) ;
    LabelFont=PlotInfo.LabelFont{i,:};
    LabelBold=PlotInfo.LabelBold(i,:) ;
    if LabelBold==1
        Bold='bold';
    else
        Bold='normal';
    end
    LabelItalics=PlotInfo.LabelItalics(i,:); 
    if LabelItalics==1
        Italics='italic';
    else
        Italics='normal';
    end    
    LabelInterpreter=PlotInfo.LabelInterpreter{i,:} ;     
    for j = 1:size(AllXs,2)                 
        yPos=min(ylim)-((range(ylim)/40)*Staggers(i,j));
        text(AllXs(i,j),yPos,labels{1,j},'HorizontalAlignment','center','FontSize',LabelFontSize,'FontName',LabelFont,'FontWeight',Bold,'FontAngle',Italics,'Interpreter',LabelInterpreter,'Rotation',LabelRotate);
        hold on
    end
end
if IncludeLegend==1
    legend(UseAxes,'Location',LegendLocation)
end

numYSlots=factorial(height(PlotInfo)-1);
YUnit=range(ylim)/100;
if strcmpi(UsePs,'Ttest')
    [~,~,PairwisePs,IndvPs] = PairwiseDiffPermStats(PlotInfo.Data,10000,XLinePosition,'Two-tailed');
else
    [PairwisePs,IndvPs] = PairwiseDiffPermStats(PlotInfo.Data,10000,XLinePosition,'Two-tailed');
end
[ ~,IndvPs ] = SigSeg(IndvPs);
%Plot Indv astricies
for i = 1:size(XsForStats,2)
    for j = 1:size(XsForStats,1)
        if ~isempty(IndvPs{j,i})
            UseX=XsForStats(j,i);
            UseY=MaxPts(j,i)+YUnit*2;
            text(UseX,UseY,IndvPs{j,i},'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Interpreter','none');
            hold on    
        end
    end
end

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
                        text(UseX,UseYSig,tempPairwisePs{j,k},'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Interpreter','none');
                        hold on    
                        plot([UseX1,UseX2],[UseYLine,UseYLine],'k');
                        hold on
                    end
                end
            end
        end
    end    

end
if IncludeLegend==1
    legend(UseAxes(1:height(PlotInfo)),'Location',LegendLocation,'FontName',LabelFont)
end
hold off
end