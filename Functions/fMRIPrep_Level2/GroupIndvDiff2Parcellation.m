function [ParcelMap,Q,parcelRSMs] = GroupIndvDiff2Parcellation(ExperimentsDir,fmriprep_table,varargin)
%Written by David Rothlein
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[ThreshMap] = VariableSetter('ThreshMap',[],varargin);
[compiledData,~,~,~,~,~,fmriprep_table_name,LoadParams] = CompileND(ExperimentsDir,fmriprep_table);
GroupDir=strrep(ExperimentsDir,'Experiments/','IndvDiff_GroupAnalysis/');
GroupDir=strrep(GroupDir,'//','/');
try
    GroupDir=[GroupDir,fmriprep_table_name,'/',LoadParams.AnalysisType,'/',LoadParams.AnalysisName,'/',LoadParams.ParcelName{1,1},'/'];
catch
    GroupDir=[GroupDir,fmriprep_table_name,'/',LoadParams.AnalysisType,'/',LoadParams.AnalysisName,'/',LoadParams.ParcelName,'/'];
end
compiledData(compiledData==0)=nan;
overlapMap=sum(single(~isnan(compiledData)),4);
load('DefaultMasks/HalkoAnat_ProbSeg_3mmMNI.mat','Anat') 
nThreshold=uiEnterName('',['Enter refRSM overlap theshold out of ', num2str(max(overlapMap(:)))]);
nThreshold=str2num(nThreshold);
if ~exist(GroupDir,'file')
    mkdir(GroupDir);     
end
saveName=uiEnterName('IndvDiffParcel','Enter Save name');
pThresh=uiEnterName('','Set p threshold (leave blank for none).');
AnatMask=single(Anat.CSF<0.15 & Anat.GM > 0.05);

overlapMap=repmat(overlapMap,[1,1,1,size(compiledData,4)]);
overlapMap=overlapMap>nThreshold;
compiledData(overlapMap==0)=nan;
compiledData=reshape(compiledData,[size(compiledData,1)*size(compiledData,2)*size(compiledData,3),size(compiledData,4)]);
useMask=overlapMap(:,:,:,1);
if isempty(ThreshMap)
    ThreshMap=single(useMask);
else
    ThreshMap=single(ThreshMap);
end
ParcelMask=single(useMask).*single(AnatMask).*ThreshMap;
ParcelVec=ParcelMask(:);

useRSMs=compiledData(ParcelVec==1,:);
useRSMs=fillmissing(useRSMs,'constant',nanmean(useRSMs,2),2);
useRSMs=useRSMs';

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
export_fig([GroupDir,saveName,'_GroupedMatrix.png'],'-png','-m5','-transparent');
close

figure
histogram(groupVec)
export_fig([GroupDir,saveName,'_GroupHistogram.png'],'-png','-m5','-transparent');
close

SaveBrik_3mmMNI(ParcelMap,{'ParcelMap'},[GroupDir,saveName,'_BrainMap']);
parcelRSMs=zeros(size(useRSMs,1),max(groupVec));
for j = 1:max(groupVec)
    parcelRSMs(:,j)=nanmean(useRSMs(:,groupVec==j),2);
end
save([GroupDir,saveName],'parcelRSMs','ParcelMask','Q','groupVec','ParcelMap');
end

