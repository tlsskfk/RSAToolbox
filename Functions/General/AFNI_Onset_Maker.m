function AFNI_Onset_Maker( SaveName,Data,NewLine)
%Written by David Rothlein
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if NewLine==1
    fileID = fopen(SaveName,'a');
    fprintf(fileID,'\r\n');
else
    fileID = fopen(SaveName,'w');
end
if isempty(Data)
    fprintf(fileID,'%c','*');
elseif size(Data,1)==1
	fprintf(fileID,'%.2f\t',Data);
elseif size(Data,1)==2
    fprintf(fileID,'%.2f*%.2f\t',Data);  
elseif size(Data,1)==3    
    fprintf(fileID,'%.2f*%.2f,%.2f\t',Data);  
else
    fprintf(fileID,'%.2f*%.2f,%.2f,%.2f\t',Data);   
end
    
fclose(fileID);
