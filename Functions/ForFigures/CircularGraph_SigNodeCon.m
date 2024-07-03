function [grph,EmptyGraph] = CircularGraph_SigNodeCon(rVec,pVec,rMat,pMat,Labels,varargin)
%Written by David Rothlein to work with CircularGraph
bonf=0.05/((length(rMat)^2-length(rMat))/2);
bonf2=0.05/length(Labels);

EmptyGraph=1;
pRSM_pVec = VariableSetter( 'pRSM_pVec',[],varargin);
pRSM_rVec = VariableSetter( 'pRSM_rVec',[],varargin);
pRSM_Labels = VariableSetter( 'pRSM_Labels',[],varargin);

params.FontSize = VariableSetter( 'FontSize',14,varargin);

params.FontBold = VariableSetter( 'FontBold','normal',varargin); %'normal' or 'bold'
params.FontItalic = VariableSetter( 'FontItalic','normal',varargin); %'normal' or 'italic'
params.pThresh = VariableSetter( 'pThresh',[bonf],varargin); 
params.pThresh=[params.pThresh,bonf];
params.pThreshVec = VariableSetter( 'pThreshVec',[bonf2],varargin);
params.pThreshVec=[params.pThreshVec,bonf2];
params.binerize = VariableSetter( 'Binarize',1,varargin); 
params.scale = VariableSetter('scale',{'sig'},varargin);
params.Sign = VariableSetter( 'Sign','Both',varargin);
%params.pThresh [0.05,0.01, ...]
%params.tails   1 or 2
%params.binarize   1 or 0
%params.scale [min,max] or [] or 'sig'
sigType=[{'-.'},{'-'}];
%params.pThresh=[0.1,0.05,0.005];
numOfNodes=length(rMat);
if isempty(Labels)
    for i = 1:numOfNodes
        Labels{1,i}=num2str(i);
    end
end

if iscell(rMat)
    [ TestMat ] = CellConvert( rMat );
    [ TestSig ] = CellConvert( pMat );
else
    TestMat=rMat;
    TestSig=pMat;
end    
if iscell(rVec)
    [ rVec ] = cell2mat( rVec );
    [ pVec ] = cell2mat( pVec );
end    
pSig=pVec<=params.pThreshVec(1,1);
if any(pSig,'all')
    EmptyGraph=0;
end
vecSign=sign(rVec);
sigVec=single(pSig).*single(vecSign);
sigVec=sigVec(:);
    colorVec=cell(length(sigVec),1);
for i = 1:length(colorVec)
    if sigVec(i,1)==1
        colorVec{i,1}=[1.0000,0.6275,0.4784];
        Labels{i,1}=[Labels{i,1},'*'];
    elseif sigVec(i,1) == -1
        colorVec{i,1}=[0.2745,0.5098,0.7059];
        Labels{i,1}=[Labels{i,1},'*'];
    else
        colorVec{i,1}=[0,0,0];
    end
end
pSig=pVec<=params.pThreshVec(1,2);
if any(pSig,'all')
    EmptyGraph=0;
end
sigVec=single(pSig);
sigVec=sigVec(:);

boldVec=zeros(length(sigVec),1);
for i = 1:length(boldVec)
    if sigVec(i,1)==1
        Labels{i,1}=[Labels{i,1},'*'];
        boldVec(i,1)=1;
    end
end
italVec=zeros(length(sigVec),1);
pRSM_Append=cell(size(italVec,1),1);
if ~isempty(pRSM_pVec)
    pRSM_Labels=pRSM_Labels(:);
    if iscell(pRSM_rVec)
        [ rVec ] = cell2mat( pRSM_rVec );
        [ pVec ] = cell2mat( pRSM_pVec );
    else
        rVec = pRSM_rVec;
        pVec = pRSM_pVec;
    end
    pSig=single(pVec<=params.pThreshVec(1,1))+single(pVec<=params.pThreshVec(1,2));
    if any(pSig~=0,'all')
        EmptyGraph=0;
    end
    signVec=sign(rVec);
    for i = 1:size(rVec,2)
        try
            tempLabel=pRSM_Labels{i,1}(1,1:3);
        catch
            tempLabel=pRSM_Labels{i,1};
        end
        for j = 1:size(rVec,1)
            if pSig(j,i)==1
                if signVec(j,i)==1
                    pRSM_Append{j,1}=[pRSM_Append{j,1},'\fontsize{10} \bf \color{orange} ^{',tempLabel,'*} '];
                else
                    pRSM_Append{j,1}=[pRSM_Append{j,1},'\fontsize{10} \bf \color{lightBlue} ^{',tempLabel,'*} '];
                end
            elseif pSig(j,i)==2
                if signVec(j,i)==1
                    pRSM_Append{j,1}=[pRSM_Append{j,1},'\fontsize{10} \bf \color{orange} ^{',tempLabel,'**} '];
                else
                    pRSM_Append{j,1}=[pRSM_Append{j,1},'\fontsize{10} \bf \color{lightBlue} ^{',tempLabel,'**} '];
                end
            end
        end
    end
end
[ TestSig2 ] = SigSeg2( TestSig,params.pThresh,'Pos');
TestSig=TestSig2~=0;
if any(TestSig,'all')
    EmptyGraph=0;
end
Labels=Labels(:);


TestMatPos=TestMat;
TestMatPos(TestMat<0)=nan;
TestMatPos(isnan(TestMatPos))=0;
TestMatPos=TestMatPos.*single(TestSig);
if params.binerize == 1
    TestMatPos(TestMatPos~=0)=1;
end

TestMatNeg=TestMat;
TestMatNeg(TestMat>=0)=nan;
TestMatNeg(isnan(TestMatNeg))=0;
TestMatNeg=TestMatNeg.*single(TestSig);
if params.binerize == 1
    TestMatNeg(TestMatNeg~=0)=1;
end 

TestMatZero=TestMatNeg*0;
run=1;
if strcmpi(params.Sign,'both')
    if sum(TestMatPos(:))==0
        if sum(TestMatNeg(:))==0
            run=0;
        else
            params.Sign='neg';
        end
    elseif sum(TestMatNeg(:))==0
        params.Sign='pos';
    end
elseif strcmpi(params.Sign,'pos') && sum(TestMatPos(:))==0
    run=0;
elseif strcmpi(params.Sign,'neg') && sum(TestMatNeg(:))==0
    run=0;
end

if run==1
    if strcmpi(params.Sign,'both')
        BothVals=abs(TestMat(TestSig==1));
        if isempty(params.scale)
            if max(BothVals)>1
                [ BothVals ] = scaleVals( BothVals,min(BothVals),1 );
            end
            if min(BothVals)<0
                [ BothVals ] = scaleVals( BothVals,0,max(BothVals));
            end    
            TestMatBothColor=TestMat.*single(TestSig);
            TestMatBothColor(TestSig)=BothVals;
        elseif iscell(params.scale)
            TestMatBothColor=cell(1);
            TestMatBothType=cell(1);
            for n = 1:length(TestSig)
                for p = 1:length(TestSig)
                    if TestMatPos(n,p)==1 
                       %TestMatBothColor{n,p}=params.scale{1,TestSig2(n,p)};
                        TestMatBothColor{n,p}=[1.0000,0.6275,0.4784];
                        TestMatBothType{n,p}=sigType{1,TestSig2(n,p)};
                    end
                    if TestMatNeg(n,p)==1 
                        TestMatBothColor{n,p}=[0.2745,0.5098,0.7059];
                        TestMatBothType{n,p}=sigType{1,TestSig2(n,p)};
                    end                    
                end
            end            
        else
            [ BothVals ] = scaleVals( BothVals,params.scale(1,1),params.scale(1,2));  
            TestMatBothColor=TestMat.*single(TestSig);
            TestMatBothColor(TestSig)=BothVals; 
        end 
        grph=circularGraph(TestMatPos,'Label',repmat({''},numOfNodes,1),'ColorMap',repmat([.8,.8,0],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
        for n=1:length(grph.Node)
            if length(grph.Node(1,n).Connection)>1
                count=2;
                for p = 1:length(grph.Node)
                    if TestMatPos(n,p)==1
                        lineType=TestMatBothType{n,p};
                        if iscell(params.scale)
                            rgbVal=TestMatBothColor{n,p};
                        else
                            rgbVal=[1,TestMatBothColor(n,p),0.2];
                        end
                        grph.Node(1,n).Connection(1,count).LineStyle=lineType;
                        grph.Node(1,n).Connection(1,count).Color=rgbVal; 
                        count=count+1;
                    end
                end
            end
        end
        hold on
        grph=circularGraph(TestMatNeg,'Label',repmat({''},numOfNodes,1),'ColorMap',repmat([0,0,.5],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
        for n=1:length(grph.Node)
            if length(grph.Node(1,n).Connection)>1
                count=1;
                for p = 1:length(grph.Node)
                    if TestMatNeg(n,p)==1
                        lineType=TestMatBothType{n,p};
                        if iscell(params.scale)
                            rgbVal=TestMatBothColor{n,p};
                        else
                            rgbVal=[1,TestMatBothColor(n,p),0.2];
                        end
                        grph.Node(1,n).Connection(1,count).LineStyle=lineType;
                        grph.Node(1,n).Connection(1,count).Color=rgbVal; 
                        count=count+1;
                    end
                end
            end
        end        
        hold on
        grph=circularGraph(TestMatZero,'Label',Labels,'ColorMap',repmat([0,0,0],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
    elseif strcmpi(params.Sign,'pos')
        PosVals=TestMat(TestMatPos==1);
        if isempty(params.scale)==1
            if max(PosVals)>1
                [ PosVals ] = scaleVals( PosVals,min(PosVals),1 );
            end
            if min(PosVals)<0
                [ PosVals ] = scaleVals( PosVals,0,max(PosVals));
            end    
            TestMatPosColor=TestMat.*TestMatPos;
            TestMatPosColor(TestMatPos==1)=PosVals;
        elseif iscell(params.scale)
            TestMatPosColor=cell(1);
            TestMatPosType=cell(1);
            for n = 1:length(TestSig)
                for p = 1:length(TestSig)
                    if TestMatPos(n,p)==1
                        TestMatPosColor{n,p}=[1.0000,0.6275,0.4784];
                        TestMatPosType{n,p}=sigType{1,TestSig2(n,p)};
                    end
                end
            end    
        else
            [ PosVals ] = scaleVals( PosVals,params.scale(1,1),params.scale(1,2));  
            TestMatPosColor=TestMat.*TestMatPos;
            TestMatPosColor(TestMatPos==1)=PosVals;   
        end 
        grph=circularGraph(TestMatPos,'Label',repmat({''},numOfNodes,1),'ColorMap',repmat([.8,.8,0],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
        for n=1:length(grph.Node)
            if length(grph.Node(1,n).Connection)>1
                count=2;
                for p = 1:length(grph.Node)
                    if TestMatPos(n,p)==1
                        lineType=TestMatPosType{n,p};
                        if iscell(params.scale)
                            rgbVal=TestMatPosColor{n,p};
                        else                        
                            rgbVal=[1,TestMatPosColor(n,p),0.2];
                        end
                        grph.Node(1,n).Connection(1,count).LineStyle=lineType;
                        grph.Node(1,n).Connection(1,count).Color=rgbVal; 
                        count=count+1;
                    end
                end
            end
        end        
        hold on  
        grph=circularGraph(TestMatZero,'Label',Labels,'ColorMap',repmat([0,0,0],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
    elseif strcmpi(params.Sign,'neg')
        NegVals=abs(TestMat(TestMatNeg==1));
        if isempty(params.scale)==1
            if max(NegVals)>1
                [ NegVals ] = scaleVals( NegVals,min(NegVals),1 );
            end
            if min(NegVals)<0
                [ NegVals ] = scaleVals( NegVals,0,max(NegVals));
            end    
            TestMatNegColor=TestMat.*TestMatNeg;
            TestMatNegColor(TestMatNeg==1)=NegVals;
        elseif iscell(params.scale)
            TestMatNegColor=cell(1);
            TestMatNegType=cell(1);
            for n = 1:length(TestSig)
                for p = 1:length(TestSig)
                    if TestMatNeg(n,p)==1
                        TestMatNegColor{n,p}=[0.2745,0.5098,0.7059];
                        TestMatNegType{n,p}=sigType{1,TestSig2(n,p)};
                    end
                end
            end   
        else
            [ NegVals ] = scaleVals( NegVals,params.scale(1,1),params.scale(1,2));  
            TestMatNegColor=TestMat.*TestMatNeg;
            TestMatNegColor(TestMatNeg==1)=NegVals;  
        end       
        grph=circularGraph(TestMatNeg,'Label',repmat({''},numOfNodes,1),'ColorMap',repmat([0,0,.5],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';
        for n=1:length(grph.Node)
            if length(grph.Node(1,n).Connection)>1
                count=2;
                for p = 1:length(grph.Node)
                    if TestMatNeg(n,p)==1
                        lineType=TestMatNegType{n,p};
                        if iscell(params.scale)
                            rgbVal=TestMatNegColor{n,p};
                        else                        
                            rgbVal=[0,TestMatNegColor(n,p),1];
                        end
                        grph.Node(1,n).Connection(1,count).LineStyle=lineType;
                        grph.Node(1,n).Connection(1,count).Color=rgbVal; 
                        count=count+1;
                    end
                end
            end
        end  
        hold on 
        grph=circularGraph(TestMatZero,'Label',Labels,'ColorMap',repmat([0,0,0],numOfNodes,1));
        grph.HideButton.Visible='off';
        grph.ShowButton.Visible='off';        
    end
    
    for i = 1:length(grph.Node(1,1).Connection.Parent.Children)
        
        if isprop(grph.Node(1,1).Connection.Parent.Children(i,1),'Text')
            txtObj=grph.Node(1,1).Connection.Parent.Children(i,1);
            for j = 1:length(Labels)
                if strcmp(Labels{j,1},txtObj.String)

                    grph.Node(1,1).Connection.Parent.Children(i,1).Color=colorVec{j,1};
                    if ~isempty(pRSM_Append{j,1}) 
                        grph.Node(1,1).Connection.Parent.Children(i,1).String=[grph.Node(1,1).Connection.Parent.Children(i,1).String,newline,pRSM_Append{j,1}];
                    end
                    if isempty(params.FontSize)==0
                        if length(params.FontSize)==1
                            grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=params.FontSize;
                        else
                            grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=params.FontSize(1,j);
                        end
                    end
                    if boldVec(j,1) == 1  %'normal' or 'bold'
                        grph.Node(1,1).Connection.Parent.Children(i,1).FontWeight='bold';
                        %grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=grph.Node(1,1).Connection.Parent.Children(i,1).FontSize+2
                    end
                    if italVec(j,1) == 1 %'normal' or 'italic'
                        grph.Node(1,1).Connection.Parent.Children(i,1).FontAngle='italic';
                    end
                    continue
                end
            end
        end
    end
else
    grph=circularGraph(TestMatZero,'Label',Labels,'ColorMap',repmat([0,0,0],numOfNodes,1));
    grph.HideButton.Visible='off';
    grph.ShowButton.Visible='off';  
    for i = 1:length(grph.Node(1,1).Connection.Parent.Children)
        
        if isprop(grph.Node(1,1).Connection.Parent.Children(i,1),'Text')
            txtObj=grph.Node(1,1).Connection.Parent.Children(i,1);
            for j = 1:length(Labels)
                if strcmp(Labels{j,1},txtObj.String)

                    grph.Node(1,1).Connection.Parent.Children(i,1).Color=colorVec{j,1};
                    if ~isempty(pRSM_Append{j,1}) 
                        grph.Node(1,1).Connection.Parent.Children(i,1).String=[grph.Node(1,1).Connection.Parent.Children(i,1).String,newline,pRSM_Append{j,1}];
                    end
            
                    if isempty(params.FontSize)==0
                        if length(params.FontSize)==1
                            grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=params.FontSize;
                        else
                            grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=params.FontSize(1,j);
                        end
                    end
                    if boldVec(j,1) == 1  %'normal' or 'bold'
                        grph.Node(1,1).Connection.Parent.Children(i,1).FontWeight='bold';
                        %grph.Node(1,1).Connection.Parent.Children(i,1).FontSize=grph.Node(1,1).Connection.Parent.Children(i,1).FontSize+2
                    end
                    if italVec(j,1) == 1 %'normal' or 'italic'
                        grph.Node(1,1).Connection.Parent.Children(i,1).FontAngle='italic';
                    end
                    continue
                end
            end
        end
    end    
end
