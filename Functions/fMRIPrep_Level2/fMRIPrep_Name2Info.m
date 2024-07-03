function [fMRIPrepInfo] = fMRIPrep_Name2Info(fMRIPrep_Name)
%Written by David Rothlein
%Convert a BIDs formatted file name into a struct containing file properties

temp_uBreaks=strfind(fMRIPrep_Name,'_'); %locate all underscores in filename
hBreaks=strfind(fMRIPrep_Name,'-'); %locate hyphens in filename

%underscores could indicate a new property or they could indicate a space
%in a single property value. This code identfies those indicating new
%properties by identify cases when multiple underscores occur between
%hyphens and removing all but the last underscore in such cases.
uBreaks=zeros(1,length(hBreaks)-1);
for i = 1:length(hBreaks)-1
    % identify underscore(s) between hyphens
    uInds = temp_uBreaks(1,temp_uBreaks>hBreaks(1,i) & temp_uBreaks<hBreaks(1,i+1));
    % select index of nearest underscore preceding a hyphen
    uBreaks(1,i)=uInds(1,end);
end
    

%locate end of filename
endBreak=strfind(fMRIPrep_Name,'.');
%if filename doesn't have a .ext, the end break is the length of filename
if isempty(endBreak)
    endBreak=length(fMRIPrep_Name);
end    
endBreak=endBreak(1,1);

%Select property names and property values
for i = 1:length(hBreaks)
    %select property name
    if i==1 
        %for first hyphen, property name is btw beginning and first hyphen
        propertyName=fMRIPrep_Name(1,1:hBreaks(1,i)-1);
    else
        %property name is btw previous underscore and current hyphen
        propertyName=fMRIPrep_Name(1,uBreaks(1,i-1)+1:hBreaks(1,i)-1);
    end
    
    %select corresponding property val
    if i == length(hBreaks)
        %for last hyphen, property val is btw last hyphen and end
        propertyVal=fMRIPrep_Name(1,hBreaks(1,i)+1:endBreak);
    else
        %property val is btw previous current hyphen and current underscore
        propertyVal=fMRIPrep_Name(1,hBreaks(1,i)+1:uBreaks(1,i)-1);
    end
    
    %Create output struct where fieldname is propertyName and field val is
    %propertyVal
    propertyVal=strrep(propertyVal,'.','');
    fMRIPrepInfo.(propertyName)=propertyVal;
end
fMRIPrepInfo=struct2table(fMRIPrepInfo);
end


