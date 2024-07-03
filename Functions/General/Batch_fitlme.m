function [ Out ] = Batch_fitlme( t,Yind,input,formulae )
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    formulae=[];
end
if isempty(formulae)
    count=1;
    for i = Yind
        disp(count/length(Yind))
        Out.models{1,count}=fitlme(t,LME_formulaGen( input.Labels,i,input.FFXind,input.RFXind,input.FFXbyRFXind,'Full',input.Full,'IntRFX',input.IntRFXind,'IntFFX',input.IntFFXind,'endAppend',input.endAppend ));
        Out.yNames{1,count}=input.Labels{i,1};
        Out.coefs(:,count)=Out.models{1,count}.Coefficients.Estimate;
        Out.coef_Ts(:,count)=Out.models{1,count}.Coefficients.tStat;
        Out.coef_Ps(:,count)=Out.models{1,count}.Coefficients.pValue;
        Out.coef_names=Out.models{1,count}.Coefficients.Name;
        Out.R2(1,count)=Out.models{1,count}.Rsquared.Ordinary;
        Out.R2(2,count)=Out.models{1,count}.Rsquared.Adjusted;
        count=count+1;
    end
end



