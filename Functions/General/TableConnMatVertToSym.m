function [SymTable,SymMat] = TableConnMatVertToSym(VertTable,VariableName,labels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SymMat=vertRSM2SymRSM(VertTable.(VariableName));
%labels=VertTable.Properties.RowNames;
SymTable=array2table(SymMat,'VariableNames',labels,'RowNames',labels);
end