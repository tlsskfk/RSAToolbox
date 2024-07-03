function [GroupMatResults,Group_Mat_Ts,Group_Mat_Means,Mat_Degree_Table,ReliabilityResults] = GroupMatParcellationSummaryFigures(matVals,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [SavePrefix] = VariableSetter('SavePrefix','',varargin);
    [CMThresholds] = VariableSetter('CMThresholds',[0.05,0.01,0.005,0.001],varargin);
    [MakeFigs] = VariableSetter('MakeFigs',1,varargin);
    [ComputeReliability] = VariableSetter('ComputeReliability',1,varargin);
    [MakeBrainMaps] = VariableSetter('MakeBrainMaps',1,varargin);
    [matName] = VariableSetter('matName',['Connectivity'],varargin);

    if istable(matVals)
        matVals=table2array(matVals);
    end    
    if ndims(matVals)==3
        MatValsVert=mat2uppertriuvectormat(matVals);
    else
        MatValsVert=matVals;
        matVals = vertRSM2SymRSM( matVals);
    end
    [labelPairs1,labelPairs2,labelPairsCell]=labels2uppertriuvectorlabels(UseLabels);
    [ Mat_ind,Mat_reverseInd ] = LabelPair2Ind( UseLabels,labelPairs1,'_2_' );
    [Group_Mat_Ts,Group_Mat_Zs,Group_Mat_Means,Group_Mat_loCI,Group_Mat_hiCI,Group_Mat_p]=getTval(matVals,3);
    GroupMatResults=cat(3,Group_Mat_Means,Group_Mat_loCI,Group_Mat_hiCI,Group_Mat_Ts,Group_Mat_Zs,Group_Mat_p);
    Group_Mat_Ts=array2table(Group_Mat_Ts,'VariableNames',UseLabels,'RowNames',UseLabels);
    Group_Mat_Means=array2table(Group_Mat_Means,'VariableNames',UseLabels,'RowNames',UseLabels);
    GroupMatResults=mat2uppertriuvectormat(GroupMatResults);
    GroupMatResults=array2table(GroupMatResults,'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','p'});
    GroupMatResults=[cell2table(labelPairsCell,'VariableNames',{'ROI1','ROI2'}),GroupMatResults,cell2table([labelPairs1,labelPairs2],'VariableNames',{'LabelPairs1','LabelPairs2'}),cell2table([ Mat_ind,Mat_reverseInd ],'VariableNames',{'MatInd','ReverseMatInd'})];
    if ComputeReliability==1
        [ ReliabilityResults.Mat_SplitHalfReliability] = CrossValCorr( MatValsVert,MatValsVert,'CrossValType','SplitHalf','Reps',100);
        [ ReliabilityResults.Mat_losoReliability] = CrossValCorr( MatValsVert,MatValsVert,'CrossValType','leave1out');   
    else
        ReliabilityResults=[];
    end
    if MakeFigs == 1
        if ComputeReliability==1
            figure
                subplot(1,2,1)
                histogram(ReliabilityResults.Mat_SplitHalfReliability)
                hold on
                title([matName,' Split Half Reliability = ',num2str(mean(ReliabilityResults.Mat_SplitHalfReliability(:)),2)],'FontSize',10)
                subplot(1,2,2)
                histogram(ReliabilityResults.Mat_losoReliability)
                hold on
                title([matName,' LOSO Reliability = ',num2str(mean(ReliabilityResults.Mat_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.Mat_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.Mat_losoReliability(:),1),3)],'FontSize',8)  
                export_fig([GroupFiguresDir,SavePrefix,'Reliability_Histograms.png'],'-Transparent','-png','-m5')
            close
        end
        Group_Vec_p=ones(size(matVals,1),1);
        Group_Vec_Means=Group_Vec_p;
        try
            if size(Group_Vec_p(:),1)<50   
                if size(Group_Vec_p(:),1) < 15
                    FontSize = 12;
                elseif size(Group_Vec_p(:),1) >= 15 && size(Group_Vec_p(:),1) < 25
                    FontSize = 10;
                else
                    FontSize = 8;
                end        
                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_p05.png'],'-Transparent','-png','-m5')
                    end
                close   

                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_p005.png'],'-Transparent','-png','-m5')
                    end    
                close

                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize);
                    if EmptyGraph==0
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_bonf.png'],'-Transparent','-png','-m5')
                    end
                close                      
            end    
        end
    end
    MapLabels=[];
    mapVals=[];    
    for pThresh=1:4
        [Results] = ConnectivityMatrixStats(table2array(Group_Mat_Means),'pmat',Group_Mat_p,'Thresholds',CMThresholds(1,pThresh));
        tempVarNames=Results.Properties.VariableNames(:)';
        for vNum = 1:length(tempVarNames)
            tempVarNames{1,vNum}=[tempVarNames{1,vNum},'_',num2str4filename(CMThresholds(1,pThresh),4)];
        end
        MapLabels=[MapLabels,tempVarNames];
        mapVals=[mapVals,table2array(Results)];
    end
    mapVals(isinf(mapVals))=0;
    mapVals(isnan(mapVals))=0;
    [ ConnectivityMaps ] = MakeParcellationMap( mapVals,UseMask);
    mapName=[matName,'_',parcelName];
    if MakeBrainMaps == 1
        SaveBrik_3mmMNI(ConnectivityMaps,MapLabels,[GroupBrainMapsDir,mapName]);
    end
    save([GroupBrainMapsDir,mapName,'.mat'],'ConnectivityMaps');            
    Mat_Degree_Table=array2table(mapVals,'VariableNames',MapLabels);
    Mat_Degree_Table.Properties.RowNames=UseLabels;  
end

