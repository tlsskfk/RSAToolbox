function [Edge_Degree_Table] = IndvDiffParcellationSummaryFigures(ResultsNodes,ResultsEdges,NodeVals,EdgeVals,BehavVals,BehLabel,NodeLabel,EdgeLabel,GroupBrainMapsDir,GroupFiguresDir,parcelName,UseMask,UseLabels,varargin)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [SavePrefix] = VariableSetter('SavePrefix','',varargin);
    [CMThresholds] = VariableSetter('CMThresholds',[0.05,0.01,0.005,0.001],varargin);
    GroupingVar = VariableSetter('GroupingVar',[],varargin);
    [OnlyLMEScatter] = VariableSetter('OnlyLMEScatter',1,varargin);
    UseLabels=UseLabels(:);
    if ~isempty(ResultsNodes)
        if istable(NodeVals)
            NodeVals=table2array(NodeVals);
        end
        NodePs_lme=ResultsNodes.lme_pValue;
        Bonf_NodePs_lme = find(NodePs_lme <= 0.05/length(NodePs_lme));
        NodeVal_lme=ResultsNodes.lme_Estimate;
        NodeZs_lme=p2z(NodePs_lme,2).*sign(NodeVal_lme);
        if ~isempty(Bonf_NodePs_lme)
            for i = Bonf_NodePs_lme(:)'
                NodeName=UseLabels{i,1};
                ScatterTitle=['LME: ',NodeLabel,newline,NodeName];
                SaveName=[NodeLabel,'_',NodeName];
                if length(SaveName)>45
                    SaveName=SaveName(1,1:45);
                end                
                statsText=['LME slope = ',num2str(NodeVal_lme(i,1),4),'; Z = ',num2str(NodeZs_lme(i,1),3),'; N = ',num2str(length(BehavVals))];
                x = NodeVals(:,i);
                lmem_scatter(GroupingVar,BehavVals,x,BehLabel,NodeLabel,ScatterTitle,statsText);
                SaveDir=[GroupFiguresDir,'LME_',BehLabel,'/'];
                if ~exist(SaveDir)
                    mkdir(SaveDir)
                end      
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                close
            end
        end
        if OnlyLMEScatter==0
            NodePs_r=ResultsNodes.pearson_p;
            Bonf_NodePs_r = find(NodePs_r <= 0.05/length(NodePs_r));        
            NodeVal_r=ResultsNodes.pearson_r;
            NodeZs_r=p2z(NodePs_r,2).*sign(NodeVal_r);
            if ~isempty(Bonf_NodePs_r)
                for i = Bonf_NodePs_r(:)'
                    NodeName=UseLabels{i,1};
                    ScatterTitle=['Pearson r: ',NodeLabel,newline,NodeName];
                    SaveName=[NodeLabel,'_',NodeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['r = ',num2str(NodeVal_r(i,1),4),'; Z = ',num2str(NodeZs_r(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = NodeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,NodeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Pearson_r_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end   
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end

            NodePs_rho=ResultsNodes.spearman_p;
            Bonf_NodePs_rho = find(NodePs_rho <= 0.05/length(NodePs_rho));        
            NodeVal_rho=ResultsNodes.spearman_r;
            NodeZs_rho=p2z(NodePs_rho,2).*sign(NodeVal_rho);       
            if ~isempty(Bonf_NodePs_rho)
                for i = Bonf_NodePs_rho(:)'
                    NodeName=UseLabels{i,1};
                    ScatterTitle=['Spearman rho: ',NodeLabel,newline,NodeName];
                    SaveName=[NodeLabel,'_',NodeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['Spearman rho = ',num2str(NodeVal_rho(i,1),4),'; Z = ',num2str(NodeZs_rho(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = NodeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,NodeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Spearman_rho_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end     
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end

            NodePs_reg=ResultsNodes.reg_t_pval;
            Bonf_NodePs_reg = find(NodePs_reg <= 0.05/length(NodePs_reg));        
            NodeVal_reg=ResultsNodes.reg_beta;
            NodeZs_reg=p2z(NodePs_reg,2).*sign(NodeVal_reg);          
            if ~isempty(Bonf_NodePs_reg)
                for i = Bonf_NodePs_reg(:)'
                    NodeName=UseLabels{i,1};
                    ScatterTitle=['Reg Coef: ',NodeLabel,newline,NodeName];
                    SaveName=[NodeLabel,'_',NodeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['Reg coef = ',num2str(NodeVal_reg(i,1),4),'; Z = ',num2str(NodeZs_reg(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = NodeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,NodeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Reg_Coef_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end      
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end
        end
    end
    if ~isempty(ResultsEdges)
        labelPairs1=labels2uppertriuvectorlabels(UseLabels);
        labelPairs1=labelPairs1(:);
        if istable(EdgeVals)
            EdgeVals=table2array(EdgeVals);
        end        
        EdgePs_lme=vertRSM2SymRSM(ResultsEdges.lme_pValue);
        Bonf_EdgePs_lme = find(ResultsEdges.lme_pValue <= 0.05/length(ResultsEdges.lme_pValue));
        EdgeVal_lme=vertRSM2SymRSM(ResultsEdges.lme_Estimate);
        EdgeZs_lme=p2z(EdgePs_lme,2).*sign(EdgeVal_lme);
        if ~isempty(Bonf_EdgePs_lme)
            tempEdgeZs_lme=mat2uppertriuvectormat(EdgeZs_lme);
            for i = Bonf_EdgePs_lme(:)'
                EdgeName=labelPairs1{i,1};
                ScatterTitle=['LME: ',EdgeLabel,newline,EdgeName];
                SaveName=[EdgeLabel,'_',EdgeName];
                if length(SaveName)>45
                    SaveName=SaveName(1,1:45);
                end                
                statsText=['LME slope = ',num2str(ResultsEdges.lme_Estimate(i,1),4),'; Z = ',num2str(tempEdgeZs_lme(i,1),3),'; N = ',num2str(length(BehavVals))];
                x = EdgeVals(:,i);
                lmem_scatter(GroupingVar,BehavVals,x,BehLabel,EdgeLabel,ScatterTitle,statsText);
                SaveDir=[GroupFiguresDir,'LME_',BehLabel,'/'];
                if ~exist(SaveDir)
                    mkdir(SaveDir)
                end      
                set(gcf, 'Position', get(0, 'Screensize'));
                export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                close
            end
        end
        if OnlyLMEScatter==0
            EdgePs_r=vertRSM2SymRSM(ResultsEdges.pearson_p);
            Bonf_EdgePs_r = find(ResultsEdges.pearson_p <= 0.05/length(ResultsEdges.pearson_p)); 
            EdgeVal_r=vertRSM2SymRSM(ResultsEdges.pearson_r);
            EdgeZs_r=p2z(EdgePs_r,2).*sign(EdgeVal_r);
            if ~isempty(Bonf_EdgePs_r)
                tempEdgeZs_r=mat2uppertriuvectormat(EdgeZs_r);
                for i = Bonf_EdgePs_r(:)'
                    EdgeName=labelPairs1{i,1};
                    ScatterTitle=['Pearson r: ',EdgeLabel,newline,EdgeName];
                    SaveName=[EdgeLabel,'_',EdgeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['r = ',num2str(ResultsEdges.pearson_r(i,1),4),'; Z = ',num2str(tempEdgeZs_r(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = EdgeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,EdgeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Pearson_r_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end       
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end

            EdgePs_rho=vertRSM2SymRSM(ResultsEdges.spearman_p);
            Bonf_EdgePs_rho = find(ResultsEdges.spearman_p <= 0.05/length(ResultsEdges.spearman_p)); 
            EdgeVal_rho=vertRSM2SymRSM(ResultsEdges.spearman_r);
            EdgeZs_rho=p2z(EdgePs_rho,2).*sign(EdgeVal_rho);       
            if ~isempty(Bonf_EdgePs_rho)
                tempEdgeZs_rho=mat2uppertriuvectormat(EdgeZs_rho);
                for i = Bonf_EdgePs_rho(:)'
                    EdgeName=labelPairs1{i,1};
                    ScatterTitle=['Spearman rho: ',EdgeLabel,newline,EdgeName];
                    SaveName=[EdgeLabel,'_',EdgeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['Spearman rho = ',num2str(ResultsEdges.spearman_r(i,1),4),'; Z = ',num2str(tempEdgeZs_rho(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = EdgeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,EdgeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Spearman_rho_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end   
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end

            EdgePs_reg=vertRSM2SymRSM(ResultsEdges.reg_t_pval);
            Bonf_EdgePs_reg = find(ResultsEdges.reg_t_pval <= 0.05/length(ResultsEdges.reg_t_pval));   
            EdgeVal_reg=vertRSM2SymRSM(ResultsEdges.reg_beta);
            EdgeZs_reg=p2z(EdgePs_reg,2).*sign(EdgeVal_reg); 
            if ~isempty(Bonf_EdgePs_reg)
                tempEdgeZs_reg=mat2uppertriuvectormat(EdgeZs_reg);
                for i = Bonf_EdgePs_reg(:)'
                    EdgeName=labelPairs1{i,1};
                    ScatterTitle=['Reg coef: ',EdgeLabel,newline,EdgeName];
                    SaveName=[EdgeLabel,'_',EdgeName];
                    if length(SaveName)>45
                        SaveName=SaveName(1,1:45);
                    end                
                    statsText=['Reg coef = ',num2str(ResultsEdges.reg_beta(i,1),4),'; Z = ',num2str(tempEdgeZs_reg(i,1),3),'; N = ',num2str(length(BehavVals))];
                    x = EdgeVals(:,i);
                    lmem_scatter(GroupingVar,BehavVals,x,BehLabel,EdgeLabel,ScatterTitle,statsText);
                    SaveDir=[GroupFiguresDir,'Reg_coef_',BehLabel,'/'];
                    if ~exist(SaveDir)
                        mkdir(SaveDir)
                    end     
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SaveName,'.png'],'-Transparent','-png','-m2')
                    close
                end
            end 
        end    
        if size(UseLabels(:),1)<50   
            if size(UseLabels(:),1) < 15
                FontSize = 12;
            elseif size(UseLabels(:),1) >= 15 && size(UseLabels(:),1) < 25
                FontSize = 10;
            else
                FontSize = 8;
            end  
            if isempty(ResultsNodes)
                Group_Vec_Means=ones(length(UseLabels),1);
                Group_Vec_p=Group_Vec_Means;
            else
                Group_Vec_Means=NodeVal_lme;
                Group_Vec_p=NodePs_lme;
            end
            Group_Mat_Means=EdgeVal_lme;
            Group_Mat_p=EdgePs_lme;
            SaveDir=[GroupFiguresDir,'LME_',BehLabel,'/'];
                if ~exist(SaveDir)
                    mkdir(SaveDir)
                end   
            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.05,'pThreshVec',0.05);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SavePrefix,'CircPlot_p05.png'],'-Transparent','-png','-m2')
                end
            close   

            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize,'pThresh',0.005,'pThreshVec',0.005);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SavePrefix,'CircPlot_p005.png'],'-Transparent','-png','-m2')
                end    
            close

            figure
                [~,EmptyGraph]=CircularGraph_SigNodeCon(Group_Vec_Means(:),Group_Vec_p(:),Group_Mat_Means,Group_Mat_p,strrepCell(UseLabels,'_',' '),'FontSize',FontSize);
                if EmptyGraph==0
                    set(gcf, 'Position', get(0, 'Screensize'));
                    export_fig([SaveDir,SavePrefix,'CircPlot_bonf.png'],'-Transparent','-png','-m2')
                end
            close     
        end         
    end    
     
  
    if ~isempty(ResultsNodes)
        SaveDir=[GroupBrainMapsDir,'LME_',BehLabel,'/'];
        if ~exist(SaveDir)
            mkdir(SaveDir)
        end           
        Group_Vec_Means=NodeVal_lme;
        Group_Vec_Zs=NodeZs_lme;        
        mapName=[NodeLabel,'_',parcelName];
        mapVals=[Group_Vec_Means(:),Group_Vec_Zs(:)];
        mapVals(isinf(mapVals))=0;
        mapVals(isnan(mapVals))=0;
        [ ActivationMaps ] = MakeParcellationMap( mapVals,UseMask);
        SaveBrik_3mmMNI(ActivationMaps,{['Mean_',NodeLabel],['Z_',NodeLabel]},[SaveDir,mapName]);
        save([SaveDir,mapName,'.mat'],'ActivationMaps');    
    end    
    if ~isempty(ResultsEdges)
        MapLabels=[];
        mapVals=[];    
        Group_Mat_Means=EdgeVal_lme;
        Group_Mat_p=EdgePs_lme;  
        for pThresh=1:4
            try
                [Results] = ConnectivityMatrixStats(Group_Mat_Means,'pmat',Group_Mat_p,'Thresholds',CMThresholds(1,pThresh));
            catch
                disp('Error')
                continue
            end
            tempVarNames=Results.Properties.VariableNames(:)';
            for vNum = 1:length(tempVarNames)
                tempVarNames{1,vNum}=[tempVarNames{1,vNum},'_',num2str4filename(CMThresholds(1,pThresh),4)];
            end
            MapLabels=[MapLabels,tempVarNames];
            mapVals=[mapVals,table2array(Results)];
        end
        SaveDir=[GroupBrainMapsDir,'LME_',BehLabel,'/'];
        if ~exist(SaveDir)
            mkdir(SaveDir)
        end                 
        mapVals(isinf(mapVals))=0;
        mapVals(isnan(mapVals))=0;
        [ ConnectivityMaps ] = MakeParcellationMap( mapVals,UseMask);
        mapName=[EdgeLabel,'_',parcelName];
        SaveBrik_3mmMNI(ConnectivityMaps,MapLabels,[SaveDir,mapName]);
        save([SaveDir,mapName,'.mat'],'ConnectivityMaps');            
        Edge_Degree_Table=array2table(mapVals,'VariableNames',MapLabels);
        Edge_Degree_Table.Properties.RowNames=UseLabels;  
    else
        Edge_Degree_Table=[];
    end
end

