function [ Behavior ] = gradCPT_SplitHalfBehavior( SplitTCbySS,SplitTCbySS_smooth,inFilesBySS,computeFull,useLabels,computeSplit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    useLabels=[];
    computeFull=0;
    computeSplit=1;
end
if nargin==4
    useLabels=[];
    computeSplit=1;
end
if nargin==5
    computeSplit=1;
end

load('gradCPT/ind_var_scenes_mix');
ind_var_city=ind_var_city(:,[1,6,11]);
ind_var_mount=ind_var_mount(:,[1,6,11]);
ind_var_labels=ind_var_labels([1,6,11]',:);
load('gradCPT/gradCPT_Output_Labels');

VTCbySS=inFilesBySS.VTC;
Zone=inFilesBySS.Zone;
ZonePrime=inFilesBySS.ZonePrime;
VTCderiv=inFilesBySS.VTCderiv;
ComputeReward=inFilesBySS.ComputeReward;



if computeFull == 1
    [Full] = GradCPT_bySubj( inFilesBySS,ind_var_city,ind_var_mount,[],VTCbySS);    
    Output=[];
    CityRT_byItem=[];
    MountainCE_byItem=[];
    MountainN_byItem=[];
    

        CityItems=Full.CityItems;
        CityRTs=Full.CityRTs;
        CityRT_byItem=Full.CityRT_byItem;
        MountainCE=Full.MountainCE;
        MountainItems=Full.MountainItems;
        Output=Full.Output;
        MountainCE_byItem=Full.MountainCE_byItem;
        MountainN_byItem=Full.MountainN_byItem;

    Output=Output(:,[6,8,7,3]);
    GroupCERate=sum(MountainCE_byItem,2)./sum(MountainN_byItem,2);
    [ Regress.Reg_RTz_byItem_Grp,Regress.Reg_CEz_byItem_Grp,Regress.Reg_RTz_byItem_SS,Regress.Reg_RTz_byTrial_SS,Regress.Reg_CEz_byTrial_SS,Regress.IndvDiff_Reg_CEz_byTrial,Regress.IndvDiff_Reg_RTz_byTrial,Regress.IndvDiff_Corr_CEz_byTrial,Regress.IndvDiff_Corr_RTz_byTrial,Regress.RTz_Index_b,RTz_Index_r,Regress.RTz_Index_t, Regress.CEz_Index_b,Regress.CEz_Index_t,Regress.CEz_Index_r ,]=GradCPT_RegressAnalysis(CityRT_byItem,{CityRTs},GroupCERate,{MountainCE},ind_var_city,ind_var_mount,{CityItems},{MountainItems},Output,'linear');
    Full.Output_Full=[Full.Output_Full,Regress.Reg_RTz_byTrial_SS.t,Regress.Reg_RTz_byTrial_SS.R,Regress.Reg_CEz_byTrial_SS.t,Regress.Reg_CEz_byTrial_SS.R];
    Behavior.FullResults=Full;
    
    if ~isempty(useLabels)
        Behavior.FullOutputLabels=useLabels;
        Behavior.FullOutput=nan(1,length(useLabels));
        if size(useLabels,1)<length(useLabels)
            useLabels=useLabels';
        end
        for n=1:length(useLabels)
            iUse=[];
            for lbln=1:length(gradCPT_Output_Labels)
                if strcmp(useLabels{n,1},gradCPT_Output_Labels{1,lbln})
                    iUse=lbln;
                    break
                end
            end
            Behavior.FullOutput(1,n)=Full.Output_Full(1,iUse);

        end
    else
        Behavior.FullOutput=nan(1,length(gradCPT_Output_Labels));
        Behavior.FullOutputLabels=gradCPT_Output_Labels;
        nlength=min([length(gradCPT_Output_Labels);length(Full.Output_Full)]);
        for n=1:nlength
            Behavior.FullOutput(1,n)=Full.Output_Full(1,n);  
        end
    end
end

if computeSplit==1
    [Split1] = GradCPT_bySubj( inFilesBySS,ind_var_city,ind_var_mount,SplitTCbySS_smooth==1,VTCbySS);
    [Split2] = GradCPT_bySubj( inFilesBySS,ind_var_city,ind_var_mount,SplitTCbySS_smooth==0,VTCbySS);



    Output=[];
    CityRT_byItem=[];
    MountainCE_byItem=[];
    MountainN_byItem=[];
    CityItems=Split1.CityItems;
    CityRTs=Split1.CityRTs;
    CityRT_byItem=Split1.CityRT_byItem;
    MountainCE=Split1.MountainCE;
    MountainItems=Split1.MountainItems;
    Output=Split1.Output;
    MountainCE_byItem=Split1.MountainCE_byItem;
    MountainN_byItem=Split1.MountainN_byItem; 
    Output=Output(:,[6,8,7,3]);
    GroupCERate=sum(MountainCE_byItem,2)./sum(MountainN_byItem,2);
    [ Regress.Reg_RTz_byItem_Grp,Regress.Reg_CEz_byItem_Grp,Regress.Reg_RTz_byItem_SS,Regress.Reg_RTz_byTrial_SS,Regress.Reg_CEz_byTrial_SS,Regress.IndvDiff_Reg_CEz_byTrial,Regress.IndvDiff_Reg_RTz_byTrial,Regress.IndvDiff_Corr_CEz_byTrial,Regress.IndvDiff_Corr_RTz_byTrial,Regress.RTz_Index_b,RTz_Index_r,Regress.RTz_Index_t, Regress.CEz_Index_b,Regress.CEz_Index_t,Regress.CEz_Index_r]=GradCPT_RegressAnalysis(CityRT_byItem,{CityRTs},GroupCERate,{MountainCE},ind_var_city,ind_var_mount,{CityItems},{MountainItems},Output,'linear');
    Split1.Output_Full=[Split1.Output_Full,Regress.Reg_RTz_byTrial_SS.t,Regress.Reg_RTz_byTrial_SS.R,Regress.Reg_CEz_byTrial_SS.t,Regress.Reg_CEz_byTrial_SS.R];


    Output=[];
    CityRT_byItem=[];
    MountainCE_byItem=[];
    MountainN_byItem=[];
    CityItems=Split2.CityItems;
    CityRTs=Split2.CityRTs;
    CityRT_byItem=Split2.CityRT_byItem;
    MountainCE=Split2.MountainCE;
    MountainItems=Split2.MountainItems;
    Output=Split2.Output;
    MountainCE_byItem=Split2.MountainCE_byItem;
    MountainN_byItem=Split2.MountainN_byItem;
    Output=Output(:,[6,8,7,3]);
    GroupCERate=sum(MountainCE_byItem,2)./sum(MountainN_byItem,2);
    [ Regress.Reg_RTz_byItem_Grp,Regress.Reg_CEz_byItem_Grp,Regress.Reg_RTz_byItem_SS,Regress.Reg_RTz_byTrial_SS,Regress.Reg_CEz_byTrial_SS,Regress.IndvDiff_Reg_CEz_byTrial,Regress.IndvDiff_Reg_RTz_byTrial,Regress.IndvDiff_Corr_CEz_byTrial,Regress.IndvDiff_Corr_RTz_byTrial,Regress.RTz_Index_b,RTz_Index_r,Regress.RTz_Index_t, Regress.CEz_Index_b,Regress.CEz_Index_t,Regress.CEz_Index_r]=GradCPT_RegressAnalysis(CityRT_byItem,{CityRTs},GroupCERate,{MountainCE},ind_var_city,ind_var_mount,{CityItems},{MountainItems},Output,'linear');
    Split2.Output_Full=[Split2.Output_Full,Regress.Reg_RTz_byTrial_SS.t,Regress.Reg_RTz_byTrial_SS.R,Regress.Reg_CEz_byTrial_SS.t,Regress.Reg_CEz_byTrial_SS.R];


    if ~isempty(useLabels)
        if size(useLabels,1)<length(useLabels)
            useLabels=useLabels';
        end
        for n=1:length(useLabels)
            iUse=[];
            for lbln=1:length(gradCPT_Output_Labels)
                if strcmp(useLabels{n,1},gradCPT_Output_Labels{1,lbln})
                    iUse=lbln;
                    break
                end
            end

            SplitMats.Split1.(gradCPT_Output_Labels{1,iUse})=Split1.Output_Full(1,iUse);
            SplitMats.Split2.(gradCPT_Output_Labels{1,iUse})=Split2.Output_Full(1,iUse);
            SplitMats.SplitDiff.(gradCPT_Output_Labels{1,iUse})=Split1.Output_Full(1,iUse)-Split2.Output_Full(1,iUse);

        end
    else
        for n=1:length(gradCPT_Output_Labels)
            SplitMats.Split1.(gradCPT_Output_Labels{1,n})=Split1.Output_Full(1,n);
            SplitMats.Split2.(gradCPT_Output_Labels{1,n})=Split2.Output_Full(1,n);
            SplitMats.SplitDiff.(gradCPT_Output_Labels{1,n})=Split1.Output_Full(1,n)-Split2.Output_Full(1,n);  
        end
    end

    TC=SplitTCbySS;
    Ind1=TC==1;
    Ind2=TC==0;
    VTC=VTCbySS;
    ZoneDeriv=ZonePrime;
    if ComputeReward ==1
        RewardTC=single(inFilesBySS.bordertracker(:,2)==255);
    end
    N=sum(single(TC~=0));
    SplitMats.Split1.meanVTC=nanmean(VTC(Ind1,1));
    SplitMats.Split2.meanVTC=nanmean(VTC(Ind2,1));
    SplitMats.SplitDiff.meanVTC=SplitMats.Split1.meanVTC-SplitMats.Split2.meanVTC;

    SplitMats.Split1.medianVTC=nanmedian(VTC(Ind1,1));
    SplitMats.Split2.medianVTC=nanmedian(VTC(Ind2,1));
    SplitMats.SplitDiff.medianVTC=SplitMats.Split1.medianVTC-SplitMats.Split2.medianVTC;

    SplitMats.Split1.meanVTCderiv=nanmean(VTCderiv(Ind1,1));
    SplitMats.Split2.meanVTCderiv=nanmean(VTCderiv(Ind2,1));
    SplitMats.SplitDiff.meanVTCderiv=SplitMats.Split1.meanVTCderiv-SplitMats.Split2.meanVTCderiv;

    SplitMats.Split1.medianVTCderiv=nanmedian(VTCderiv(Ind1,1));
    SplitMats.Split2.medianVTCderiv=nanmedian(VTCderiv(Ind2,1));
    SplitMats.SplitDiff.medianVTCderiv=SplitMats.Split1.medianVTCderiv-SplitMats.Split2.medianVTCderiv;

    SplitMats.Split1.ZonePercent=sum(Zone(Ind1,1))/sum(Ind1);
    SplitMats.Split2.ZonePercent=sum(Zone(Ind2,1))/sum(Ind2);
    SplitMats.SplitDiff.ZonePercent=SplitMats.Split1.ZonePercent-SplitMats.Split2.ZonePercent;

    SplitMats.Split1.ZoneDerivPercent=sum(ZoneDeriv(Ind1,1))/sum(Ind1);
    SplitMats.Split2.ZoneDerivPercent=sum(ZoneDeriv(Ind2,1))/sum(Ind2);
    SplitMats.SplitDiff.ZoneDerivPercent=SplitMats.Split1.ZoneDerivPercent-SplitMats.Split2.ZoneDerivPercent;
    if ComputeReward ==1
        SplitMats.Split1.RewardPercent=sum(RewardTC(Ind1,1))/sum(Ind1);
        SplitMats.Split2.RewardPercent=sum(RewardTC(Ind2,1))/sum(Ind2);
        SplitMats.SplitDiff.RewardPercent=SplitMats.Split1.RewardPercent-SplitMats.Split2.RewardPercent;
    end

    % 
    % AllFields=fields(SplitMats.Split1);
    % numFields=size(AllFields,1);
    % 
    % for n = 1:numFields
    %     for j = 1
    %         TCmat=zeros(numSplits,length(SplitTCbySS{1,1}));
    %         for i = 1:numSplits
    %             tempC=SplitTCbySS_smooth{i,1}';
    %             TCmat(i,tempC==1)=SplitMats.Split1.(AllFields{n,1})(i,1);
    %             TCmat(i,tempC==2)=SplitMats.Split2.(AllFields{n,1})(i,1);
    %         end
    %         TCs.(AllFields{n,1}){1,1}=nanmean(TCmat,1);
    %         TC_Zs.(AllFields{n,1}){1,1}=zscore(nanmean(TCmat,1));
    %     end
    % end
    % TCs.origVTC=VTC;
    % TC_Zs.origVTC=zscore(VTC);
    % TCs.origVTCderiv=VTCderiv;
    % TC_Zs.origVTCderiv=zscore(VTCderiv);
    % TCs.origZone=Zone;
    % TCs.origZoneDeriv=ZoneDeriv;
    % TC_Zs.origZone=zscore(Zone);
    % TC_Zs.origZoneDeriv=zscore(ZoneDeriv);
    % 
    % if ComputeReward ==1
    %     TCs.RewardTC=RewardTC;
    %     TC_Zs.RewardTC=zscore(RewardTC); 
    % end


    Behavior.SplitMats=SplitMats;
    % Behavior.TCs=TCs;
    % Behavior.TC_Zs=TC_Zs;

    Behavior.SplitTCbySS=SplitTCbySS;
    Behavior.SplitTCbySS_smooth=SplitTCbySS_smooth;
end
