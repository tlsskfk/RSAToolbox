function [ParcelMap,Q,parcelRSMs] = GroupRCSearchlight2Parcellation(ExperimentsDir,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ThreshMap] = VariableSetter('ThreshMap',[],varargin);

[~,~,fmriprep_table_name] = load_fmriprep_table('ExperimentsDir',ExperimentsDir);
startDir=strrep(ExperimentsDir,'Experiments',['GroupAnalysis/',fmriprep_table_name,'/RSMs/']);
useDir=uigetdir(startDir);
load([useDir,'/','slRSMs.mat'],'OverlapVec','slMask','slRSMs','slZRSMs');
load('DefaultMasks/HalkoAnat_ProbSeg_3mmMNI.mat','Anat') 
nThreshold=uiEnterName('',['Enter refRSM overlap theshold out of ', num2str(max(OverlapVec))]);
nThreshold=str2num(nThreshold);
if isempty(ThreshMap)
    ThreshMap=single(slMask);
else
    ThreshMap=single(ThreshMap);
end
saveName=uiEnterName('','Enter Save name');
pThresh=uiEnterName('','Set p threshold (leave blank for none).');
AnatMask=single(Anat.CSF<0.15 & Anat.GM > 0.05);
AnatVec=AnatMask(slMask);
ThreshVec=ThreshMap(slMask);
OverlapMap=single(slMask)*0;
OverlapMap(slMask)=OverlapVec;
OverlapMap=single(OverlapMap >= nThreshold);
ParcelMask=single(AnatMask).*single(OverlapMap).*single(ThreshMap);
useVec=single(OverlapVec >= nThreshold).*single(AnatVec).*single(ThreshVec);
useRSMs=slRSMs(:,useVec==1);
%slZRSMs=slZRSMs(:,useVec==1);
if ~isempty(pThresh)
    pThresh=str2num(pThresh);
    [r,p]=corrcoef(useRSMs);
    r(p<pThresh)=0;
else
    r=corrcoef(useRSMs);
end

[groupVec,Q]=community_louvain(r,[],[],'negative_asym');
ParcelMap=ParcelMask;
ParcelMap(ParcelMask==1)=groupVec;
[~,i]=sort(groupVec);
CommunityMat=r(i,i);
figure
imagesc(CommunityMat)
export_fig([useDir,'/',saveName,'_GroupedMatrix.png'],'-png','-m5','-transparent');
close

figure
histogram(groupVec)
export_fig([useDir,'/',saveName,'_GroupHistogram.png'],'-png','-m5','-transparent');
close

SaveBrik_3mmMNI(ParcelMap,{'ParcelMap'},[useDir,'/',saveName,'_BrainMap.png']);
parcelRSMs=zeros(size(useRSMs,1),max(groupVec));
for j = 1:max(groupVec)
    parcelRSMs(:,j)=nanmean(useRSMs(:,groupVec==j),2);
end
save([useDir,'/',saveName],'parcelRSMs','ParcelMask','Q','groupVec','ParcelMap');
end

