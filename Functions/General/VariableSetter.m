function [ UseVar,varName ] = VariableSetter( varName,defaultVal,varArgNames1,varArgVals2)
%Written by David Rothlein
% functions called: vararginConvert


if nargin==3
    [ varArgNames1,varArgVals2 ] = vararginConvert( varArgNames1 );
end
numVarArgs=length(varArgNames1);
varFind=0;
for i = 1:numVarArgs
    varArgName=varArgNames1{1,i};
    varVal=varArgVals2{1,i};
    if strcmpi(varName,varArgName)==1
        UseVar=varVal;
        varFind=1;
        break
    end   
end

if varFind==0
   UseVar = defaultVal;
end
end

function [ VarOutName,VarOutVal ] = vararginConvert( VarCell )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numArg=size(VarCell,2);
VarOutName=cell(1,4);
VarOutVal=cell(1,4);
countA=1;
countB=1;
for i= 1:numArg
    if mod(i,2)==1
        VarOutName{1,countA}=VarCell{1,i};
        countA=countA+1;
    else
        VarOutVal{1,countB}=VarCell{1,i};
        countB=countB+1;
    end
end

end