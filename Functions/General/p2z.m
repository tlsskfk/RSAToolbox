function [zvals]=p2z(pvals,tails)
%Written by David Rothlein
if tails==1

    zvals=icdf('normal',(1-(pvals)),0,1);    
    
else
    
    zvals=icdf('normal',(1-(pvals/2)),0,1);

end
   