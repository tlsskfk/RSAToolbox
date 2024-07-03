function SmoothedData=bvsmoother(data,fwhm) 

    %data is 3D matrix
    %fwhm = full width half mass (gaussian decay rate)

    sigma_val=fwhm2sigma(fwhm);
    if fwhm<1
        kernelsize=1;
    else
        kernelsize=(fwhm*2)+1;
    end

    SmoothedData=smooth3(data,'gaussian',[kernelsize,kernelsize,kernelsize],sigma_val);
end

function fwhm_sigma = fwhm2sigma(fwhm)

    fwhm_sigma = fwhm/(sqrt(8*log(2)));
end