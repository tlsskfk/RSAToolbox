function [RSMs,ImgCoords] = ImageRSMSearchlight(ImageMat,WindowSize)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
numImg=size(ImageMat,3);
ImgH=size(ImageMat,1);
ImgW=size(ImageMat,2);
NumSLs=(ImgH-WindowSize)*(ImgW-WindowSize);
RSMs=nan((numImg^2-numImg)/2,NumSLs);
ImgCoords=nan(NumSLs,2);
count=1;
for i = 1:(ImgH-WindowSize)
    for j = 1:(ImgW-WindowSize)
        %disp(count/NumSLs)
       ImgCoords(count,:)=[i,j];
       ImgSeg=ImageMat(i:i+WindowSize-1,j:j+WindowSize-1,:);
       ImgSeg=reshape(ImgSeg,[WindowSize^2,numImg]);
       filt=sum(isnan(ImgSeg),2)==numImg;
       if sum(single(filt))/length(filt) < 0.5
            RSMs(:,count)=mat2uppertriuvectormat(corr(ImgSeg,'rows','pairwise'));
       end
       count=count+1;
    end
end
ImgCoords=ImgCoords+round(WindowSize/2);
end

