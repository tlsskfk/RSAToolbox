function img=ImageSimMatImg(RSM,cmapType,range,DataTxt,ImgPaths,imScaleX,imScaleY )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    range=[];
    DataTxt=0;
    ImgPaths=[];
    cmapType='gray';
end

if nargin==2
    range=[];
    DataTxt=0;
    ImgPaths=[];
end

if nargin==3
    DataTxt=0;
    ImgPaths=[];
end

if nargin==4
    ImgPaths=[];
end

rsmSize=length(RSM);
fig1=figure; 
set(fig1,'Position',[300,0,800,600],'Resize','off');

matsize=size(RSM,1);          %# A 5-by-5 matrix of random values from 0 to 1
triuTemp=triu(ones(rsmSize));
RSM(triuTemp==1)=nan;
if ~isempty(range)
    RSM(1,1)=range(1,1); 
    RSM(2,2)=range(1,2); 
end
colormap(cmapType);  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
imagesc(RSM,'AlphaData',single(triuTemp==0));            %# Create a colored plot of the matrix values
if ~isempty(range)
    RSM(1,1)=nan; 
    RSM(2,2)=nan; 
end
if DataTxt==1                         
    textStrings = num2str(RSM(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    for i = 1:length(textStrings)
        if strcmp(textStrings{i,1},'NaN')
            textStrings{i,1}='';
        end
    end

    [x,y] = meshgrid(1:matsize);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center','FontUnits','points','FontSize',12);
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(RSM(:) < midValue,1,3);  %# Choose white or black for the
                                                 %#   text color of the strings so
                                                 %#   they can be easily seen over
                                                 %#   the background color
    %%%                                             
    textColors = textColors*0+1;     
    %%%
    set(hStrings,{'Color'},num2cell(textColors,2));                          
end
a=get(gca,'position') ;
b=xlim;
c=ylim;
axis off
if ~isempty(ImgPaths)
    for j=1:rsmSize
        img=imread(ImgPaths{j,1});
        imgcoordsX=plotimage(a,0,rsmSize-j+1,b,c,[imScaleX, imScaleY]);
        imgcoordsY=plotimage(a,j,0,b,c,[imScaleX, imScaleY]);
        axes('pos',[imgcoordsX(1,1) imgcoordsX(1,2) imScaleX imScaleY]);
        ppp=imshow(img);
        mask=img~=255;
        set(ppp,'AlphaData',single(mask));
        axis off
        axes('pos',[imgcoordsY(1,1) imgcoordsY(1,2) imScaleX imScaleY]);
        ppp=imshow(img);
        set(ppp,'AlphaData',single(mask));
        axis off
    end
end

F = getframe(gcf);
[img] = frame2im(F);
close

