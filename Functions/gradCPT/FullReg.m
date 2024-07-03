function [ Stats ] = FullReg( yVars,xVars,yLabels,xLabels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sumLabs={'betas';'norm_betas';'ts';'IndvRsqr';'RsqrContrib'};

yVarsZ=yVars;
for i = 1:size(yVars,2)
    yVarsZ(:,i)=zscore(yVars(:,i));
end;    

xVarsZ=xVars;
for i = 1:size(xVars,2)
    xVarsZ(:,i)=zscore(xVars(:,i));
end;      
    
Stats=struct();
Stats.XCorrMat=cell(size(xVars,2)+1);
Stats.XCorrMat(2:end,1)=xLabels';
Stats.XCorrMat(1,2:end)=xLabels;
Stats.XCorrMat(2:end,2:end)=num2cell(corrcoef(xVars));

Stats.YCorrMat=cell(size(yVars,2)+1);
Stats.YCorrMat(2:end,1)=yLabels';
Stats.YCorrMat(1,2:end)=yLabels;
Stats.YCorrMat(2:end,2:end)=num2cell(corrcoef(yVars));

Stats.VIF=cell(2,size(xVars,2));
Stats.VIF(1,:)=xLabels;
for i = 1:size(yVars,2)
    i*10
    yVar=yVars(:,i);
    yVarZ=yVarsZ(:,i);
    
    unnorm=regstats(yVar,xVars);
    norm =regstats(yVarZ,xVarsZ);
    Stats.(yLabels{1,i}).Rsqr=unnorm.rsquare;
    Stats.(yLabels{1,i}).AdjRsqr=unnorm.adjrsquare;
    Stats.(yLabels{1,i}).Resids=unnorm.standres;
    Stats.(yLabels{1,i}).F=unnorm.fstat;
    Stats.(yLabels{1,i}).Summary=cell(6,size(xVars,2)+1);
    Stats.(yLabels{1,i}).Summary{1,1}=yLabels{1,i};
    Stats.(yLabels{1,i}).Summary(2:end,1)=sumLabs;
    Stats.(yLabels{1,i}).Summary(1,2:end)=xLabels;
    Stats.(yLabels{1,i}).Summary(2,2:end)=num2cell(unnorm.beta(2:end,1))';
    Stats.(yLabels{1,i}).Summary(3,2:end)=num2cell(norm.beta(2:end,1))';
    Stats.(yLabels{1,i}).Summary(4,2:end)=num2cell(unnorm.tstat.t(2:end,1))';
    
    
    for j = 1:size(xVars,2)
        j
        xVar=xVars(:,j);
        xNull=xVars;
        xNull(:,j)=[];
        treg1=regstats(yVar,xVar);
        treg2=regstats(yVar,xNull);
        Stats.(yLabels{1,i}).Summary{5,j+1}=treg1.rsquare;
        Stats.(yLabels{1,i}).Summary{6,j+1}=Stats.(yLabels{1,i}).Rsqr-treg2.rsquare;
        treg3=regstats(xVar,xNull);
        Stats.VIF{2,j}=1/(1-treg3.rsquare);
    end;
    
end;


