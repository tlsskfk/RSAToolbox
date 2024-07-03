function [UseCode] = StabDiff_Img2Img_CodeGen(inPath,outPath,prompt,varargin)
%PARAMs:
%txt2Img
%CFG 2-20 where 2 = mostly AI guided image and 20 = mostly prompt guided image
% --ckpt sd-v1-4.ckpt
[ddim_steps] = VariableSetter('ddim_steps',[50],varargin);
% use switch to use plms 0 1
[plms] = VariableSetter('plms',[0],varargin);
% eta=0.0 corresponds to deterministic sampling
[ddim_eta] = VariableSetter('ddim_eta',0,varargin);
% sample this often
[n_iter] = VariableSetter('n_iter',[2],varargin);
% latent channels
[C] = VariableSetter('C',4,varargin);
% downsampling factor, most often 8 or 16
[f] = VariableSetter('f',8,varargin);
% how many samples to produce for each given prompt. A.k.a batch size
[n_samples] = VariableSetter('n_samples',[3],varargin);
% unconditional guidance scale: eps = eps(x, empty) + scale * (eps(x, cond) - eps(x, empty))
[scale] = VariableSetter('scale',5,varargin);
% strength for noising/unnoising. 1.0 corresponds to full destruction of information in init image
[strength] = VariableSetter('strength',[0.75],varargin);
% the starting seed (for reproducible sampling)
[seed] = VariableSetter('seed',42,varargin);
% evaluate at this precision choices=["full", "autocast"]
[precision] = VariableSetter('precision','autocast',varargin);


UseCode = ['python scripts/img2img.py --prompt "',prompt,'" --init-img ',inPath,' --outdir ',outPath,...
    ' --ddim_steps ',num2str(ddim_steps),...
    ' --ddim_eta ',num2str(ddim_eta),...
    ' --n_iter ',num2str(n_iter),...
    ' --C ',num2str(C),...
    ' --f ',num2str(f),...
    ' --n_samples ',num2str(n_samples),...
    ' --scale ',num2str(scale),...
    ' --strength ',num2str(strength),...
    ' --seed ',num2str(seed),...
    ' --precision ',precision];

UseCode=[UseCode,' --skip_grid --ckpt sd-v1-4.ckpt '];
if plms == 1
    UseCode=[UseCode,'--plms'];
end

end

%python scripts/img2img.py --prompt "Photorealistic image of a picnic table" 
% --init-img /mnt/c/Users/david/OneDrive/matlabtoolboxes/Experiments/B_gradCPT_X5_digits_/StabDiff/orig/1.jpg 
% --ckpt sd-v1-4.ckpt --skip_grid --n_samples 3 --strength 0.95 
% --outdir /mnt/c/Users/david/OneDrive/matlabtoolboxes/Experiments/B_gradCPT_X5_digits_/StabDiff/StableDiffusion/1';