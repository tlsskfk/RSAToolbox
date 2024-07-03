function [Ts,Zs,Ms,loCI,hiCI,ps]=getTval(mat,dim)
%Written by David Rothlein
if nargin==1
    if iscell(mat)==1
        dim = length(size(mat{1,1}));
    else
        dim = length(size(mat));
    end
end
if iscell(mat)==1
    tot=length(mat);
    if size(mat,2)<tot
        mat=mat';
    end
    Ts=cell(1,tot);
    Zs=cell(1,tot);
    Ms=cell(1,tot);
    for i = 1:tot
        Ms{1,i}=nanmean(mat{1,i},dim);
        [~,ps,CI,a]=ttest(mat{1,i},0,0.05,'both',dim);        
        Ts{1,i}=a.tstat;
        Zs{1,i}=p2z(ps,2).*single(sign(Ts{1,i}));  
        loCI{1,i}=CI(1,:);
        hiCI{1,i}=CI(2,:);
    end
elseif dim == 1
    Ms=nanmean(mat,dim);
    [~,ps,CI,a]=ttest(mat,0,0.05,'both',dim);
    Ts=a.tstat;
    Zs=p2z(ps,2).*single(sign(Ts));  
    loCI=CI(1,:);
    hiCI=CI(2,:);
elseif dim == 2
    Ms=nanmean(mat,dim);
    [~,ps,CI,a]=ttest(mat,0,0.05,'both',dim);
    Ts=a.tstat;
    Zs=p2z(ps,2).*single(sign(Ts));  
    loCI=CI(:,1);
    hiCI=CI(:,2);
elseif dim == 3
    Ms=nanmean(mat,dim);
    [~,ps,CI,a]=ttest(mat,0,0.05,'both',dim);
    Ts=a.tstat;
    Zs=p2z(ps,2).*single(sign(Ts));  
    loCI=CI(:,:,1);
    hiCI=CI(:,:,2);    
elseif dim == 4
    Ms=nanmean(mat,dim);
    [~,ps,CI,a]=ttest(mat,0,0.05,'both',dim);
    Ts=a.tstat;
    Zs=p2z(ps,2).*single(sign(Ts));  
    loCI=CI(:,:,:,1);
    hiCI=CI(:,:,:,2);   
elseif dim == 5
    Ms=nanmean(mat,dim);
    [~,ps,CI,a]=ttest(mat,0,0.05,'both',dim);
    Ts=a.tstat;
    Zs=p2z(ps,2).*single(sign(Ts));  
    loCI=CI(:,:,:,:,1);
    hiCI=CI(:,:,:,:,2);      
end




    