function SymRSMs=triRSMs2SymRSMs(RSMs)
%Written by David Rothlein
SymRSMs=RSMs;

for i = 1:size(RSMs,3)
    RSM=RSMs(:,:,i);
    SymRSMs(:,:,i)=RSM+permute(RSM,[2,1]);
end
    
