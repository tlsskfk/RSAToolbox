function [ShiftedTC] = TCShift(TC,shiftVal)
%Written by David Rothlein
%Positive shiftVal = moves timepoints earlier (to the left)
%Negative shiftVal moves timepoints later (to the right)
isHorz=0;
if size(TC,2) == length(TC)
    isHorz=1;
    TC=TC(:);
end
nanpad=nan(abs(shiftVal),1);
if shiftVal > 0
    ShiftedTC=[TC(shiftVal+1:end);nanpad];
else
    ShiftedTC=[nanpad;TC(1:end+shiftVal)];
end
if isHorz==1
    ShiftedTC=ShiftedTC';
end
    


end

