function [GroupVecResults,GroupMatResults,Group_Mat_Ts,Group_Mat_Means,Mat_Degree_Table,ReliabilityResults] = GroupParcellationSummaryFigures(vecVals,matVals,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,varargin)
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
    [permZMat] = VariableSetter('permZMat',[],varargin);
    [permPMat] = VariableSetter('permPMat',[],varargin);
    [permZVec] = VariableSetter('permZVec',[],varargin);
    [permPVec] = VariableSetter('permPVec',[],varargin);
    [permZMat_Corrected] = VariableSetter('permZMat_Corrected',[],varargin);
    [permPMat_Corrected] = VariableSetter('permPMat_Corrected',[],varargin);
    [permZVec_Corrected] = VariableSetter('permZVec_Corrected',[],varargin);
    [permPVec_Corrected] = VariableSetter('permPVec_Corrected',[],varargin); 
    ReliabilityResults=[];
    if istable(vecVals)
        vecVals=table2array(vecVals);
    end
    if size(vecVals,1)~=length(UseLabels)
        vecVals=vecVals';
    end    
    if istable(matVals)
        matVals=table2array(matVals);
    end    
    if ndims(matVals)==3
        MatValsVert=mat2uppertriuvectormat(matVals);
    else
        MatValsVert=matVals;
        matVals = vertRSM2SymRSM( matVals);
    end
    [Group_Vec_Ts,Group_Vec_Zs,Group_Vec_Means,Group_Vec_loCI,Group_Vec_hiCI,Group_Vec_p]=getTval(vecVals,2); 
    if ~isempty(permZVec) && ~isempty(permZVec_Corrected)
        GroupVecResults=array2table([Group_Vec_Means(:),Group_Vec_loCI(:),Group_Vec_hiCI(:),Group_Vec_Ts(:),Group_Vec_Zs(:),permZVec(:),permZVec_Corrected(:),Group_Vec_p(:),permPVec(:),permPVec_Corrected(:)],'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Perm_Z','Corrected_Perm_Z','p','Perm_p','Corrected_Perm_p'});
    elseif ~isempty(permZVec)
        GroupVecResults=array2table([Group_Vec_Means(:),Group_Vec_loCI(:),Group_Vec_hiCI(:),Group_Vec_Ts(:),Group_Vec_Zs(:),permZVec(:),Group_Vec_p(:),permPVec(:)],'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Perm_Z','p','Perm_p'});
    elseif ~isempty(permZVec_Corrected)
        GroupVecResults=array2table([Group_Vec_Means(:),Group_Vec_loCI(:),Group_Vec_hiCI(:),Group_Vec_Ts(:),Group_Vec_Zs(:),permZVec_Corrected(:),Group_Vec_p(:),permPVec_Corrected(:)],'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Corrected_Perm_Z','p','Corrected_Perm_p'});   
    else
        GroupVecResults=array2table([Group_Vec_Means(:),Group_Vec_loCI(:),Group_Vec_hiCI(:),Group_Vec_Ts(:),Group_Vec_Zs(:),Group_Vec_p(:)],'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','p'});
    end
    GroupVecResults=[cell2table(UseLabels(:),'VariableNames',{'ROIs'}),GroupVecResults];
    [labelPairs1,labelPairs2,labelPairsCell]=labels2uppertriuvectorlabels(UseLabels);
    [ Mat_ind,Mat_reverseInd ] = LabelPair2Ind( UseLabels,labelPairs1,'_2_' );
    [Group_Mat_Ts,Group_Mat_Zs,Group_Mat_Means,Group_Mat_loCI,Group_Mat_hiCI,Group_Mat_p]=getTval(matVals,3);
    Group_Mat_Ts=array2table(Group_Mat_Ts,'VariableNames',UseLabels,'RowNames',UseLabels);
    Group_Mat_Means=array2table(Group_Mat_Means,'VariableNames',UseLabels,'RowNames',UseLabels);
    
    if ~isempty(permZMat) && ~isempty(permZMat_Corrected)
        GroupMatResults=cat(3,table2array(Group_Mat_Means),Group_Mat_loCI,Group_Mat_hiCI,table2array(Group_Mat_Ts),Group_Mat_Zs,permZMat,permZMat_Corrected,Group_Mat_p,permPMat,permPMat_Corrected);
        GroupMatResults=mat2uppertriuvectormat(GroupMatResults);
        GroupMatResults=array2table(GroupMatResults,'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Perm_Z','Corrected_Perm_Z','p','Perm_p','Corrected_Perm_p'});
    elseif ~isempty(permZMat)
        GroupMatResults=cat(3,table2array(Group_Mat_Means),Group_Mat_loCI,Group_Mat_hiCI,table2array(Group_Mat_Ts),Group_Mat_Zs,permZMat_Corrected,Group_Mat_p,permPMat_Corrected);
        GroupMatResults=mat2uppertriuvectormat(GroupMatResults);
        GroupMatResults=array2table(GroupMatResults,'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Corrected_Perm_Z','p','Corrected_Perm_p'});
    elseif ~isempty(permZMat_Corrected)
        GroupMatResults=cat(3,table2array(Group_Mat_Means),Group_Mat_loCI,Group_Mat_hiCI,table2array(Group_Mat_Ts),Group_Mat_Zs,permZMat,Group_Mat_p,permPMat);
        GroupMatResults=mat2uppertriuvectormat(GroupMatResults);
        GroupMatResults=array2table(GroupMatResults,'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','Perm_Z','p','Perm_p'});  
    else
        GroupMatResults=cat(3,table2array(Group_Mat_Means),Group_Mat_loCI,Group_Mat_hiCI,table2array(Group_Mat_Ts),Group_Mat_Zs,Group_Mat_p);
        GroupMatResults=mat2uppertriuvectormat(GroupMatResults);        
        GroupMatResults=array2table(GroupMatResults,'VariableNames',{'Mean','CI_lo','CI_hi','T','Z','p'});
    end    
    GroupMatResults=[cell2table(labelPairsCell,'VariableNames',{'ROI1','ROI2'}),GroupMatResults,cell2table([labelPairs1,labelPairs2],'VariableNames',{'LabelPairs1','LabelPairs2'}),cell2table([ Mat_ind,Mat_reverseInd ],'VariableNames',{'MatInd','ReverseMatInd'})];
    if ComputeReliability==1
        [ ReliabilityResults.Vec_SplitHalfReliability] = CrossValCorr( vecVals,vecVals,'CrossValType','SplitHalf','Reps',100);
        [ ReliabilityResults.Mat_SplitHalfReliability] = CrossValCorr( MatValsVert,MatValsVert,'CrossValType','SplitHalf','Reps',100);
        [ ReliabilityResults.Vec_losoReliability] = CrossValCorr( vecVals,vecVals,'CrossValType','leave1out');
        [ ReliabilityResults.Mat_losoReliability] = CrossValCorr( MatValsVert,MatValsVert,'CrossValType','leave1out');  
    end
    if MakeFigs==1
        if ComputeReliability==1
            figure
                subplot(2,2,1)
                histogram(ReliabilityResults.Vec_SplitHalfReliability)
                hold on
                title([vecName,' Split Half Reliability = ',num2str(mean(ReliabilityResults.Vec_SplitHalfReliability(:)),2)],'FontSize',10);
                subplot(2,2,2)
                histogram(ReliabilityResults.Mat_SplitHalfReliability)
                hold on
                title([matName,' Split Half Reliability = ',num2str(mean(ReliabilityResults.Mat_SplitHalfReliability(:)),2)],'FontSize',10)
                subplot(2,2,3)
                histogram(ReliabilityResults.Vec_losoReliability)
                hold on
                title([vecName,' LOSO Reliability = ',num2str(mean(ReliabilityResults.Vec_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.Vec_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.Vec_losoReliability(:),1),3)],'FontSize',8)
                subplot(2,2,4)
                histogram(ReliabilityResults.Mat_losoReliability)
                hold on
                title([matName,' LOSO Reliability = ',num2str(mean(ReliabilityResults.Mat_losoReliability(:)),2),'; t(',num2str(length(ReliabilityResults.Mat_losoReliability)-1),') = ',num2str(getTval(ReliabilityResults.Mat_losoReliability(:),1),3)],'FontSize',8)  
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,'Reliability_Histograms.png'],'-Transparent','-png','-m2')
            close
        end
        if length(Group_Vec_Ts)<60  
            if length(Group_Vec_Ts) < 14
                FontSize=14;
            elseif length(Group_Vec_Ts) < 24
                FontSize=11;
            elseif length(Group_Vec_Ts) < 30
                FontSize=7; 
            else
                FontSize=5;
            end
            tempVecVals=vecVals';
            PlotData=cell(1,size(tempVecVals,2));
            for i = 1:size(tempVecVals,2)
                PlotData{1,i}=tempVecVals(:,i);
            end
            labels=UseLabels(:)';
            figure
                DataPointDropPlot( PlotData,labels,'LabelFontSize',FontSize)
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([GroupFiguresDir,SavePrefix,'DataDropPlot.png'],'-Transparent','-png','-m2')
            close
        end    
        if size(Group_Vec_p(:),1)<50   
            if size(Group_Vec_p(:),1) < 15
                FontSize = 16;
            elseif size(Group_Vec_p(:),1) >= 15 && size(Group_Vec_p(:),1) < 25
                FontSize = 12;
            else
                FontSize = 10;
            end        
            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([GroupFiguresDir,SavePrefix,'CircPlot_p05.png'],'-Transparent','-png','-m2')
                end
            close   

            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([GroupFiguresDir,SavePrefix,'CircPlot_p005.png'],'-Transparent','-png','-m2')
                end    
            close

            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),table2array(Group_Mat_Means),Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([GroupFiguresDir,SavePrefix,'CircPlot_bonf.png'],'-Transparent','-png','-m2')
                end
            close  
            if ~isempty(permPMat)
                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),permPVec(:),table2array(Group_Mat_Means),permPMat,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
                    if EmptyGraph==0
                        set(gcf, 'Position', get(0, 'Screensize'));
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_perm05.png'],'-Transparent','-png','-m2')
                    end
                close   

                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),permPVec(:),table2array(Group_Mat_Means),permPMat,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
                    if EmptyGraph==0
                        set(gcf, 'Position', get(0, 'Screensize'));
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_perm005.png'],'-Transparent','-png','-m2')
                    end    
                close 
            end
            if ~isempty(permPMat_Corrected)
                figure
                    [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),permPVec_Corrected(:),table2array(Group_Mat_Means),permPMat_Corrected,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
                    if EmptyGraph==0
                        set(gcf, 'Position', get(0, 'Screensize'));
                        export_fig([GroupFiguresDir,SavePrefix,'CircPlot_perm05_Corrected.png'],'-Transparent','-png','-m2')
                    end
                close   
            end
        end
    end
    if MakeBrainMaps == 1
        if ~isempty(permZVec) && ~isempty(permZVec_Corrected)
            mapName=[vecName,'_',parcelName];
            mapVals=[Group_Vec_Means(:),Group_Vec_Ts(:),permZVec(:),permZVec_Corrected(:)];
            mapVals(isinf(mapVals))=0;
            mapVals(isnan(mapVals))=0;
            [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ActivationMaps,{['Mean_',vecName],['T_',vecName],['PermZ_',vecName],['PermZ_Corrected_',vecName]},[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps');         
        elseif ~isempty(permZVec)
            mapName=[vecName,'_',parcelName];
            mapVals=[Group_Vec_Means(:),Group_Vec_Ts(:),permZVec(:)];
            mapVals(isinf(mapVals))=0;
            mapVals(isnan(mapVals))=0;
            [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ActivationMaps,{['Mean_',vecName],['T_',vecName],['PermZ_',vecName]},[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps'); 
        elseif ~isempty(permZVec_Corrected)
            mapName=[vecName,'_',parcelName];
            mapVals=[Group_Vec_Means(:),Group_Vec_Ts(:),permZVec_Corrected(:)];
            mapVals(isinf(mapVals))=0;
            mapVals(isnan(mapVals))=0;
            [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ActivationMaps,{['Mean_',vecName],['T_',vecName],['PermZ_Corrected_',vecName]},[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps');   
        else
            mapName=[vecName,'_',parcelName];
            mapVals=[Group_Vec_Means(:),Group_Vec_Ts(:)];
            mapVals(isinf(mapVals))=0;
            mapVals(isnan(mapVals))=0;
            [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
            SaveBrik_3mmMNI(ActivationMaps,{['Mean_',vecName],['T_',vecName]},[GroupBrainMapsDir,mapName]);
            save([GroupBrainMapsDir,mapName,'.mat'],'ActivationMaps');    
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
    if ~isempty(permPMat)
        for pThresh=1:4
            [Results] = ConnectivityMatrixStats(table2array(Group_Mat_Means),'pmat',permPMat,'Thresholds',CMThresholds(1,pThresh));
            tempVarNames=Results.Properties.VariableNames(:)';
            for vNum = 1:length(tempVarNames)
                tempVarNames{1,vNum}=[tempVarNames{1,vNum},'_Perm',num2str4filename(CMThresholds(1,pThresh),4)];
            end
            MapLabels=[MapLabels,tempVarNames];
            mapVals=[mapVals,table2array(Results)];
        end 
    end
    if ~isempty(permPMat_Corrected)
        [Results] = ConnectivityMatrixStats(table2array(Group_Mat_Means),'pmat',permPMat_Corrected,'Thresholds',0.05);
        tempVarNames=Results.Properties.VariableNames(:)';
        for vNum = 1:length(tempVarNames)
            tempVarNames{1,vNum}=[tempVarNames{1,vNum},'_PermCorrected',num2str4filename(CMThresholds(1,pThresh),4)];
        end
        MapLabels=[MapLabels,tempVarNames];
        mapVals=[mapVals,table2array(Results)];        
    end
    if MakeBrainMaps == 1
        mapVals(isinf(mapVals))=0;
        mapVals(isnan(mapVals))=0;
        [ ConnectivityMaps ] = MakeParcellationMap( mapVals,UseMask);
        mapName=[matName,'_',parcelName];
        SaveBrik_3mmMNI(ConnectivityMaps,MapLabels,[GroupBrainMapsDir,mapName]);
        save([GroupBrainMapsDir,mapName,'.mat'],'ConnectivityMaps');   
    end
    Mat_Degree_Table=array2table(mapVals,'VariableNames',MapLabels);
    Mat_Degree_Table.Properties.RowNames=UseLabels;  
end

