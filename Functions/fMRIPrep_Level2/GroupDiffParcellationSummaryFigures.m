function [Mat_Degree_Table,ReliabilityResults] = GroupDiffParcellationSummaryFigures(inData_Node,inData_Edge,NodeStats,EdgeStats,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,GroupLabels,groupPairLabels,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[SavePrefix] = VariableSetter('SavePrefix','',varargin);
[MakeFigs] = VariableSetter('MakeFigs',1,varargin);
[ComputeReliability] = VariableSetter('ComputeReliability',1,varargin);
[MakeBrainMaps] = VariableSetter('MakeBrainMaps',1,varargin);
[CMThresholds] = VariableSetter('CMThresholds',[0.05,0.01,0.005,0.001],varargin);
[vecName] = VariableSetter('vecName',['Activation'],varargin);
[matName] = VariableSetter('matName',['Connectivity'],varargin);
Mat_Degree_Table=[];
ReliabilityResults=[];
numGroups=length(GroupLabels);
numGroupPairs=length(groupPairLabels);
numNodes=length(UseLabels);
numEdges=size(inData_Edge,2);

%% Create Circular Plots
if MakeFigs==1
if ~isempty(EdgeStats) && numNodes < 50
    if numNodes < 15
        FontSize = 12;
    elseif numNodes >= 15 && numNodes < 25
        FontSize = 10;
    else
        FontSize = 7;
    end  
    for i = 1:numGroups
        GroupLabel = GroupLabels{i,1};
        if ~isempty(NodeStats)
            Group_Vec_Means=NodeStats.Indv_val(i,:);
            Group_Vec_p=NodeStats.Indv_tP(i,:);
        else
            Group_Vec_Means=ones(1,numNodes);
            Group_Vec_p=Group_Vec_Means;
        end
        Group_Mat_Means=vertRSM2SymRSM(EdgeStats.Indv_val(i,:)');
        Group_Mat_p=vertRSM2SymRSM(EdgeStats.Indv_tP(i,:)');
        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'p < 0.05','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_p05.png'],'-Transparent','-png','-m2')
            end
        close   

        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'p < 0.005','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_p005.png'],'-Transparent','-png','-m2')
            end    
        close

        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'bonf(p) < 0.05','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_bonf.png'],'-Transparent','-png','-m2')
            end
        close  
    end
    
    for i = 1:numGroupPairs
        GroupLabel = groupPairLabels{i,1};
        if ~isempty(NodeStats)
            Group_Vec_Means=NodeStats.Pairwise_val(i,:);
            Group_Vec_p=NodeStats.Pairwise_lmeP(i,:);
        else
            Group_Vec_Means=ones(1,numNodes);
            Group_Vec_p=Group_Vec_Means;
        end        

        Group_Mat_Means=vertRSM2SymRSM(EdgeStats.Pairwise_val(i,:)');
        Group_Mat_p=vertRSM2SymRSM(EdgeStats.Pairwise_lmeP(i,:)');
        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'p < 0.05','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_p05.png'],'-Transparent','-png','-m2')
            end
        close   

        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'p < 0.005','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_p005.png'],'-Transparent','-png','-m2')
            end    
        close

        figure
            [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize);
            hold on
            ys=ylim;
            xs=xlim;            
            text(xs(1),ys(2),GroupLabel,'Interpreter','None','Fontsize',12,'HorizontalAlignment','right');
            text(xs(2),ys(1),'bonf(p) < 0.05','Interpreter','None','Fontsize',8);
            if EmptyGraph==0
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,GroupLabel,'_CircPlot_bonf.png'],'-Transparent','-png','-m2')
            end
        close  
    end    
end

%% Create Drop Plots
if ~all(isempty(inData_Node)) && numNodes<60 
    if numNodes < 10
        FontSize=11;
    elseif numNodes < 20
        FontSize=9;
    elseif numNodes < 30
        FontSize=7; 
    else
        FontSize=5;
    end
    PairwisePs=NodeStats.Pairwise_permP;
    IndvPs=NodeStats.Indv_permP;
    labels=UseLabels(:)';
    figure
        DataPointDropPlot( inData_Node,labels,'LabelFontSize',FontSize,'PairwisePs',PairwisePs,'IndvPs',IndvPs,'GroupNames',GroupLabels)
        set(gcf, 'Position', get(0, 'Screensize'));
        export_fig([GroupFiguresDir,vecName,'_DataDropPlot_permPs.png'],'-Transparent','-png','-m2')
    close    
    PairwisePs=NodeStats.Pairwise_lmeP;
    IndvPs=NodeStats.Indv_tP;    
    figure
        DataPointDropPlot( inData_Node,labels,'LabelFontSize',FontSize,'PairwisePs',PairwisePs,'IndvPs',IndvPs,'GroupNames',GroupLabels)
        set(gcf, 'Position', get(0, 'Screensize'));
        export_fig([GroupFiguresDir,vecName,'_DataDropPlot_lmePs.png'],'-Transparent','-png','-m2')
    close     
end    
if ~all(isempty(inData_Edge)) && numEdges<60 
    if numEdges < 10
        FontSize=11;
    elseif numEdges < 20
        FontSize=9;
    elseif numEdges < 30
        FontSize=7; 
    else
        FontSize=5;
    end
    PairwisePs=EdgeStats.Pairwise_permP;
    IndvPs=EdgeStats.Indv_permP;
    LabelPairs=labels2uppertriuvectorlabels(UseLabels);
    labels=LabelPairs(:)';
    figure
        DataPointDropPlot( inData_Edge,labels,'LabelFontSize',FontSize,'PairwisePs',PairwisePs,'IndvPs',IndvPs,'GroupNames',GroupLabels)
        set(gcf, 'Position', get(0, 'Screensize'));
        export_fig([GroupFiguresDir,matName,'DataDropPlot_permPs.png'],'-Transparent','-png','-m2')
    close    
    PairwisePs=EdgeStats.Pairwise_lmeP;
    IndvPs=EdgeStats.Indv_tP;    
    figure
        DataPointDropPlot( inData_Edge,labels,'LabelFontSize',FontSize,'PairwisePs',PairwisePs,'IndvPs',IndvPs,'GroupNames',GroupLabels)
        set(gcf, 'Position', get(0, 'Screensize'));
        export_fig([GroupFiguresDir,matName,'DataDropPlot_lmePs.png'],'-Transparent','-png','-m2')
    close     
end  
end
%% Create Reliability Plots
    if ~isempty(inData_Node)
        for i = 1:numGroups-1
            Label1=GroupLabels{i,1};
            for j = i+1:numGroups
                Label2=GroupLabels{j,1};
                LabelPair=[Label1,'vs',Label2];
                Data1=squeeze(cell2nDMAT(inData_Node(i,:)))';
                Data2=squeeze(cell2nDMAT(inData_Node(j,:)))';
                if ComputeReliability==1
                    [ ReliabilityResults.(Label1).Vec_losoReliability] = CrossValCorr( Data1,Data1,'CrossValType','leave1out');
                    [ ReliabilityResults.(Label2).Vec_losoReliability] = CrossValCorr( Data2,Data2,'CrossValType','leave1out');
                    [ ReliabilityResults.(LabelPair).Vec_losoReliability] = CrossValCorr( [Data1,Data2],[Data1,Data2],'CrossValType','leave1out');
                    if MakeFigs == 1
                        figure
                            subplot(1,3,1)
                            histogram(ReliabilityResults.(Label1).Vec_losoReliability)
                            hold on
                            title([Label1,' ',vecName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(Label1).Vec_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(Label1).Vec_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(Label1).Vec_losoReliability(:),1),3)],'FontSize',8);
                            subplot(1,3,2)
                            histogram(ReliabilityResults.(Label2).Vec_losoReliability)
                            title([Label2,' ',vecName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(Label2).Vec_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(Label2).Vec_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(Label2).Vec_losoReliability(:),1),3)],'FontSize',8);
                            subplot(1,3,3)
                            histogram(ReliabilityResults.(LabelPair).Vec_losoReliability)
                            title([LabelPair,' ',vecName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(LabelPair).Vec_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(LabelPair).Vec_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(LabelPair).Vec_losoReliability(:),1),3)],'FontSize',8);
                            set(gcf, 'Position', get(0, 'Screensize'));
                            export_fig([GroupFiguresDir,SavePrefix,LabelPair,'_Reliability_Histograms.png'],'-Transparent','-png','-m2')                  
                        close
                    end
                end
            end
        end
    end
    if ~isempty(inData_Edge)
        for i = 1:numGroups-1
            Label1=GroupLabels{i,1};
            for j = i+1:numGroups
                Label2=GroupLabels{j,1};
                LabelPair=[Label1,'vs',Label2];
                Data1=squeeze(cell2nDMAT(inData_Edge(i,:)))';
                Data2=squeeze(cell2nDMAT(inData_Edge(j,:)))';
                if ComputeReliability==1
                    [ ReliabilityResults.(Label1).Mat_losoReliability] = CrossValCorr( Data1,Data1,'CrossValType','leave1out');
                    [ ReliabilityResults.(Label2).Mat_losoReliability] = CrossValCorr( Data2,Data2,'CrossValType','leave1out');
                    [ ReliabilityResults.(LabelPair).Mat_losoReliability] = CrossValCorr( [Data1,Data2],[Data1,Data2],'CrossValType','leave1out');
                    if MakeFigs == 1
                        figure
                            subplot(1,3,1)
                            histogram(ReliabilityResults.(Label1).Mat_losoReliability)
                            hold on
                            title([Label1,' ',matName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(Label1).Mat_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(Label1).Mat_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(Label1).Mat_losoReliability(:),1),3)],'FontSize',8);
                            subplot(1,3,2)
                            histogram(ReliabilityResults.(Label2).Mat_losoReliability)
                            title([Label2,' ',matName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(Label2).Mat_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(Label2).Mat_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(Label2).Mat_losoReliability(:),1),3)],'FontSize',8);
                            subplot(1,3,3)
                            histogram(ReliabilityResults.(LabelPair).Mat_losoReliability)
                            title([LabelPair,' ',matName,' LOSO(r) = ',num2str(mean(ReliabilityResults.(LabelPair).Mat_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.(LabelPair).Mat_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.(LabelPair).Mat_losoReliability(:),1),3)],'FontSize',8);
                            set(gcf, 'Position', get(0, 'Screensize'));
                            export_fig([GroupFiguresDir,SavePrefix,LabelPair,'_Reliability_Histograms.png'],'-Transparent','-png','-m2')                    
                        close
                    end
                end
            end
        end
    end

%% Create Brain Maps 
    if ~isempty(inData_Node)
        mapName=[vecName,'_',parcelName];
        mapVals=[];
        MapLabels=[];
        for i = 1:numGroups
            Label1=GroupLabels{i,1}; 
            Group_Vec_Means=NodeStats.Indv_val(i,:)';
            Group_Vec_permZ=p2z(NodeStats.Indv_permP(i,:),2)'.*sign(Group_Vec_Means);
            Group_Vec_tZ=p2z(NodeStats.Indv_tP(i,:),2)'.*sign(Group_Vec_Means); 
            Group_Vec_tZ(isinf(Group_Vec_tZ))=NodeStats.Indv_Tval(isinf(Group_Vec_tZ))';
            mapVals=[mapVals,Group_Vec_Means,Group_Vec_permZ,Group_Vec_tZ];
            MapLabels=[MapLabels,{['Mean_',Label1],['permZ_',Label1],['tZ_',Label1]}];
        end
        for i = 1:numGroupPairs
            Label1=groupPairLabels{i,1}; 
            Group_Vec_Means=NodeStats.Pairwise_val(i,:)';
            Group_Vec_permZ=p2z(NodeStats.Pairwise_permP(i,:),2)'.*sign(Group_Vec_Means);
            Group_Vec_lmeZ=p2z(NodeStats.Pairwise_lmeP(i,:),2)'.*sign(Group_Vec_Means); 
            Group_Vec_lmeZ(isinf(Group_Vec_lmeZ))=NodeStats.Pairwise_Tval(isinf(Group_Vec_lmeZ))';
            mapVals=[mapVals,Group_Vec_Means,Group_Vec_permZ,Group_Vec_lmeZ];
            MapLabels=[MapLabels,{['Mean_',Label1],['permZ_',Label1],['lmeZ_',Label1]}];
        end     
        mapVals(isinf(mapVals))=0;
        mapVals(isnan(mapVals))=0;
        if MakeBrainMaps==1
            [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ActivationMaps,MapLabels,[GroupBrainMapsDir,mapName]);   
            save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps');
        end
    end
    
    if ~isempty(inData_Edge)
        mapName=[matName,'_',parcelName];
        mapVals=[];
        MapLabels=[];   
        for i = 1:numGroups
            Label1=GroupLabels{i,1}; 
            Group_Mat_Means=EdgeStats.Indv_val(i,:)';
            Group_Mat_p=EdgeStats.Indv_permP(i,:)';
            for pThresh=1:length(CMThresholds)
                [Results] = ConnectivityMatrixStats(Group_Mat_Means,'pmat',Group_Mat_p,'Thresholds',CMThresholds(1,pThresh));
                tempVarNames=Results.Properties.VariableNames(:)';
                for vNum = 1:length(tempVarNames)
                    tempVarNames{1,vNum}=[Label1,'_',tempVarNames{1,vNum},'_',num2str4filename(CMThresholds(1,pThresh),4)];
                end
                MapLabels=[MapLabels,tempVarNames];
                mapVals=[mapVals,table2array(Results)];
            end
        end        
        for i = 1:numGroupPairs
            Label1=groupPairLabels{i,1};
            Group_Mat_Means=EdgeStats.Pairwise_val(i,:)';
            Group_Mat_p=EdgeStats.Pairwise_lmeP(i,:)';
            for pThresh=1:length(CMThresholds)
                [Results] = ConnectivityMatrixStats(Group_Mat_Means,'pmat',Group_Mat_p,'Thresholds',CMThresholds(1,pThresh));
                tempVarNames=Results.Properties.VariableNames(:)';
                for vNum = 1:length(tempVarNames)
                    tempVarNames{1,vNum}=[Label1,'_',tempVarNames{1,vNum},'_',num2str4filename(CMThresholds(1,pThresh),4)];
                end
                MapLabels=[MapLabels,tempVarNames];
                mapVals=[mapVals,table2array(Results)];
            end
        end 

        mapVals(isinf(mapVals))=0;
        mapVals(isnan(mapVals))=0;
        if MakeBrainMaps==1
            [ ConnectivityMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ConnectivityMaps,MapLabels,[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'ConnectivityMaps');  
        end
        Mat_Degree_Table=array2table(mapVals,'VariableNames',MapLabels);
        Mat_Degree_Table.Properties.RowNames=UseLabels;  
    end
end