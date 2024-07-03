function simImg=createSimMatImage(mat,XLabels,YLabels) 


matsize=size(mat,1);          %# A 5-by-5 matrix of random values from 0 to 1
simImg=imagesc(mat);            %# Create a colored plot of the matrix values
colormap(gray);  %# Change the colormap to gray (so higher values are
caxis([-1,1]);                         %#   black and lower values are white)
colorbar;                         

textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:matsize);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center','FontUnits','points','FontSize',10);
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
% %%%                                             
% textColors = textColors*0+1;     
% %%%
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'TickLabelInterpreter','none');
set(gca,'XTick',1:matsize,...                         %# Change the axes tick marks
        'XTickLabel',XLabels,...  %#   and tick labels
        'YTick',1:matsize,...
        'YTickLabel',YLabels,...
        'TickLength',[0 0]);