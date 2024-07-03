function [GroupVecResults,ReliabilityResults] = GroupVecParcellationSummaryFigures(vecVals,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [SavePrefix] = VariableSetter('SavePrefix','',varargin);
    [vecName] = VariableSetter('vecName',['Activation'],varargin);
    if istable(vecVals)
        vecVals=table2array(vecVals);
    end
    if size(vecVals,1)~=length(UseLabels)
        vecVals=vecVals';
    end
    [Group_Vec_Ts,Group_Vec_Zs,Group_Vec_Means,Group_Vec_loCI,Group_Vec_hiCI,Group_Vec_p]=getTval(vecVals,2);    
    GroupVecResults=array2table([Group_Vec_Means(:),Group_Vec_loCI(:),Group_Vec_hiCI(:),Group_Vec_Ts(:),Group_Vec_Zs(:),Group_Vec_p(:)],'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','p'});
    GroupVecResults=[cell2table(UseLabels(:),'VariableNames',{'ROIs'}),GroupVecResults];
    [ ReliabilityResults.Vec_SplitHalfReliability] = CrossValCorr( vecVals,vecVals,'CrossValType','SplitHalf','Reps',100);
    [ ReliabilityResults.Vec_losoReliability] = CrossValCorr( vecVals,vecVals,'CrossValType','leave1out'); 
    figure
        subplot(1,2,1)
        histogram(ReliabilityResults.Vec_SplitHalfReliability)
        hold on
        title([vecName,' Split Half = ',num2str(mean(ReliabilityResults.Vec_SplitHalfReliability(:)),2)],'FontSize',10);
        subplot(1,2,2)
        histogram(ReliabilityResults.Vec_losoReliability)
        hold on
        title([vecName,' LOSO = ',num2str(mean(ReliabilityResults.Vec_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.Vec_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.Vec_losoReliability(:),1),3)],'FontSize',8)
        subtitle(SavePrefix);
        export_fig([GroupFiguresDir,SavePrefix,'Reliability_Histograms.png'],'-Transparent','-png','-m2')
    close
    
    if length(Group_Vec_Ts)<60  
        if length(Group_Vec_Ts) < 10
            FontSize=14;
        elseif length(Group_Vec_Ts) < 20
            FontSize=12;
        elseif length(Group_Vec_Ts) < 30
            FontSize=10; 
        else
            FontSize=7;
        end
        tempVecVals=vecVals';
        PlotData=cell(1,size(tempVecVals,2));
        for i = 1:size(tempVecVals,2)
            PlotData{1,i}=tempVecVals(:,i);
        end
        labels=UseLabels(:)';
        figure
            DataPointDropPlot( PlotData,labels,'LabelFontSize',FontSize)
            title(SavePrefix);
            set(gcf, 'Position', get(0, 'Screensize'));
            export_fig([GroupFiguresDir,SavePrefix,'DataDropPlot.png'],'-Transparent','-png','-m2')
        close
    end    
    
    mapName=[vecName,'_',parcelName];
    mapVals=[Group_Vec_Means(:),Group_Vec_Ts(:)];
    mapVals(isinf(mapVals))=0;
    mapVals(isnan(mapVals))=0;
    [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
    SaveBrik_3mmMNI(ActivationMaps,{['Mean_',vecName],['T_',vecName]},[GroupBrainMapsDir,mapName]);
    save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps');  
end

