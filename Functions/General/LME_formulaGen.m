function [ formula,input ] = LME_formulaGen( Labels,Yind,FFXind,RFXind,FFXbyRFXind,varargin )
%Written by David Rothlein
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[ varArgNames1,varArgVals2 ] = vararginConvert( varargin );
Full = VariableSetter( 'Full',0,varArgNames1,varArgVals2);
IntRFXind = VariableSetter( 'IntRFX',[],varArgNames1,varArgVals2);
IntFFXind = VariableSetter( 'IntFFX',[],varArgNames1,varArgVals2);
endAppend = VariableSetter( 'endAppend','',varArgNames1,varArgVals2);
input.Labels=Labels;
input.Yind=Yind;
input.FFXind=FFXind;
input.RFXind=RFXind;
input.FFXbyRFXind=FFXbyRFXind;
input.Full=Full;
input.IntRFXind=IntRFXind;
input.IntFFXind=IntFFXind;
input.endAppend=endAppend;
if size(Labels,1) ~= length(Labels)
    Labels = Labels';
end

if size(FFXind,2)~=length(FFXind)
    FFXind=FFXind';
end

Y=Labels{Yind,1};
FFX=[];
if ~isempty(FFXind)
    for i = FFXind
        FFX = [FFX,Labels{i,1},' + '];
    end
end

IntFFX=[];    
if ~isempty(IntFFXind)
    if size(IntFFXind,1)~=length(IntFFXind)
        IntFFXind=IntFFXind';
    end     
    for i = 1:size(IntFFXind,1)      
        tempInd=IntFFXind{i,1};
        temp=[];
        for j = tempInd
            if j == tempInd(1,1)
                temp=Labels{j,1};
            else
                temp=[temp,'*',Labels{j,1}];
            end
        end
        IntFFX=[IntFFX,temp,' + '];
    end
end

FFX=[FFX,IntFFX];

FFXbyRFX=[];
if Full == 1
    FFXbyRFX=FFX(1,1:end-2);
else
    if ~isempty(FFXbyRFXind)
        for i = FFXbyRFXind
            if i == FFXbyRFXind(end)
                FFXbyRFX = [FFXbyRFX,Labels{i,1}];
            else
                FFXbyRFX = [FFXbyRFX,Labels{i,1},' + '];
            end
        end
    end
    if ~isempty(IntRFXind)
        IntRFX=[];  
        if size(IntRFXind,1)~=length(IntRFXind)
            IntRFXind=IntRFXind';
        end     
        for i = 1:size(IntRFXind,1)      
            tempInd=IntRFXind{i,1};
            temp=[];
            for j = tempInd
                if j == tempInd(1,1)
                    temp=Labels{j,1};
                else
                    temp=[temp,'*',Labels{j,1}];
                end
            end
            if i == size(IntRFXind,1)
                IntRFX=[IntRFX,temp];
            else
                IntRFX=[IntRFX,temp,' + '];
            end
        end
        FFXbyRFX=[FFXbyRFX,' + ',IntRFX];
    end
end

if isempty(FFXbyRFX)
    RFX=['(1|',Labels{RFXind,1},')'];
else
    RFX=['(1 + ',FFXbyRFX,'|',Labels{RFXind,1},')'];
end

formula=[Y,' ~ ',FFX,RFX,endAppend];

