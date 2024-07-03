function [TimeSeriesParams,OutText] = uiTimeSeriesParams()
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
UIFig=uifigure('Position',[500,50,300,500],'Visible','On');
Title=uilabel(UIFig,'Position',[0,500-45,300,50],'Text','Set Timeseries Parameters','HorizontalAlignment','Center');
EndButton=uibutton(UIFig,'Position',[100,10,100,40],'Text','Done','ButtonPushedFcn',@(EndButton,event) NextButtonPush(UIFig) );

ModeTitle=uilabel(UIFig,'Position',[100,500-100,130,30],'Text',['Enter mode threshold',newline,'(between 0 and 1)'],'HorizontalAlignment','left','Visible','Off');
ModeText=uitextarea(UIFig,'Position',[20,500-100,50,30],'Value','0.5','HorizontalAlignment','Center','Visible','Off');
ModeCheck=uicheckbox(UIFig,'Position',[20,500-70,150,30],'Text','','ValueChangedFcn',@(ModeCheck,event) ModeCheckFcn(ModeTitle,ModeText) );
ModeCheckTitle=uilabel(UIFig,'Position',[100,500-70,130,30],'Text','Apply mode filter','HorizontalAlignment','left');

CVTitle=uilabel(UIFig,'Position',[100,500-170,130,30],'Text',['Enter CV threshold',newline,'(between 0 and 1)'],'HorizontalAlignment','left','Visible','Off');
CVText=uitextarea(UIFig,'Position',[20,500-170,50,30],'Value','','HorizontalAlignment','Center','Visible','Off');
CVCheck=uicheckbox(UIFig,'Position',[20,500-140,150,30],'Text','','ValueChangedFcn',@(CVCheck,event) CVCheckFcn(CVTitle,CVText) );
CVCheckTitle=uilabel(UIFig,'Position',[100,500-140,130,30],'Text','Apply CV filter','HorizontalAlignment','left');

CensorNANTitle=uilabel(UIFig,'Position',[100,500-210,200,30],'Text',['Enter value to be replaced by NAN',newline,'(leave empty to skip)'],'HorizontalAlignment','left');
CensorNANText=uitextarea(UIFig,'Position',[20,500-210,50,30],'Value','','HorizontalAlignment','Center');

ImputeTitle=uilabel(UIFig,'Position',[100,500-280,100,20],'Text','Select impute method','HorizontalAlignment','center','Visible','Off');
ImputeListBox=uilistbox(UIFig,'Position',[100,500-300,100,50],'Items',{'linear';'previous';'next';'nearest';'spline';'pchip'},'Multiselect','off','Visible','Off');
ImputeCheck=uicheckbox(UIFig,'Position',[20,500-250,200,30],'Text','','ValueChangedFcn',@(ImputeCheck,event) ImputeCheckFcn(ImputeListBox,ImputeTitle) );
ImputeCheckTitle=uilabel(UIFig,'Position',[100,500-250,130,30],'Text','Impute (fill in) NAN values','HorizontalAlignment','left');

OutlierTitle1=uilabel(UIFig,'Position',[100,500-370,150,20],'Text','Select outlier definition','HorizontalAlignment','center','Visible','Off');
OutlierListBox=uilistbox(UIFig,'Position',[100,500-390,100,40],'Items',{'mean';'median';'quartiles';'grubbs';'gesd'},'Multiselect','off','Visible','Off');
OutlierTitle2=uilabel(UIFig,'Position',[130,500-430,150,20],'Text','Enter outlier threshold','HorizontalAlignment','left','Visible','Off');
OutlierText=uitextarea(UIFig,'Position',[65,500-430,50,20],'Value','5','HorizontalAlignment','Center','Visible','Off');
OutlierCheck=uicheckbox(UIFig,'Position',[20,500-340,200,30],'Text','','ValueChangedFcn',@(OutlierCheck,event) OutlierCheckFcn(OutlierText,OutlierTitle1,OutlierTitle2,OutlierListBox) );
OutlierCheckTitle=uilabel(UIFig,'Position',[100,500-340,130,30],'Text','Detect and clip outliers','HorizontalAlignment','left');
waitfor(UIFig,'Visible','Off');

%%%% Set param values here %%%%
OutText=[];
TimeSeriesParams.Voxel_ModeRemove=single(ModeCheck.Value);
TimeSeriesParams.Voxel_ModeThresh=str2num(ModeText.Value{1,1});
if TimeSeriesParams.Voxel_ModeRemove==1
    OutText=[OutText,'Voxels were removed if the signal from that voxel remained constant for at least ',num2str(TimeSeriesParams.Voxel_ModeThresh*100,2),'% of volumes. '];
end
if strcmpi(CVText.Value{1,1},'')
    TimeSeriesParams.Voxel_CV_filter=[];
else
    TimeSeriesParams.Voxel_CV_filter=str2num(CVText.Value{1,1});
    if TimeSeriesParams.Voxel_ModeRemove==1
        OutText=[OutText,'Further, voxels were removed if the variance (CV) of their signal was greater than ', CVText.Value{1,1},'. '];
    else
        OutText=[OutText,'Voxels were removed if the variance (CV) of their signal was greater than ', CVText.Value{1,1},'. '];
    end
end

if strcmpi(CensorNANText.Value{1,1},'')
    TimeSeriesParams.Voxel_CensorNAN =[];
else
    TimeSeriesParams.Voxel_CensorNAN=str2num(CensorNANText.Value{1,1});
end
if single(ImputeCheck.Value)==0
    TimeSeriesParams.Voxel_Impute = [];
else
    TimeSeriesParams.Voxel_Impute = 1;
    if isempty(TimeSeriesParams.Voxel_CensorNAN)
        OutText=[OutText,'Within a voxel, timepoints without signal were imputed using '];
    else
        OutText=[OutText,'Within a voxel, timepoints without signal, or that had a value of ',CensorNANText.Value{1,1},' were filled in '];
    end
    if strcmpi(ImputeListBox.Value,'linear') || strcmpi(ImputeListBox.Value, 'spline')
        OutText=[OutText,'using ',ImputeListBox.Value,' interpolation. '];
    elseif strcmpi(ImputeListBox.Value,'previous') || strcmpi(ImputeListBox.Value, 'next') || strcmpi(ImputeListBox.Value, 'nearest') 
        OutText=[OutText,'by imputing the ',ImputeListBox.Value,' value with signal. '];
    else
        OutText=[OutText,'using the ',ImputeListBox.Value,' method. '];
    end
end
TimeSeriesParams.Voxel_Impute_Type = ImputeListBox.Value;

if single(OutlierCheck.Value)==0
    TimeSeriesParams.Voxel_Outlier = [];
else
    TimeSeriesParams.Voxel_Outlier = 1;
    OutText=[OutText,'Outliers were defined (relative to other timepoints within that voxel) '];
    if strcmpi(OutlierListBox.Value,'mean')
        OutText=[OutText,'as timepoints that fell outside ',OutlierText.Value{1,1},' standard deviations of the mean and were replaced with the value corresponding to the outlier threshold.'];
    elseif strcmpi(OutlierListBox.Value,'median') 
        OutText=[OutText,'as timepoints that fell outside ',OutlierText.Value{1,1},' scaled MAD from the median and were replaced with the value corresponding to the outlier threshold.'];
    else
        OutText=[OutText,'using the ',ImputeListBox.Value,' method and were replaced with the value corresponding to the outlier threshold.'];
    end 
end
TimeSeriesParams.Voxel_Outlier_Type = OutlierListBox.Value;
TimeSeriesParams.Voxel_Outlier_Thresh = str2num(OutlierText.Value{1,1});

close all force;
end

function NextButtonPush(UIFig)
    UIFig.Visible='Off';
end
function ModeCheckFcn(ModeTitle,ModeText)
    ModeTitle.Visible='On';
    ModeText.Visible='On';
end
function CVCheckFcn(CVTitle,CVText)
    CVTitle.Visible='On';
    CVText.Visible='On';
end
function ImputeCheckFcn(ImputeListBox,ImputeTitle)
    ImputeListBox.Visible='On';
    ImputeTitle.Visible='On';
end
function OutlierCheckFcn(OutlierListBox,OutlierTitle1,OutlierTitle2,OutlierText)
    OutlierListBox.Visible='On';
    OutlierTitle1.Visible='On';
    OutlierTitle2.Visible='On';
    OutlierText.Visible='On';
end
