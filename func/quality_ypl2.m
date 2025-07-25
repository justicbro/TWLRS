function [psnr, ssim, rse] = quality_ypl2(imagery1, imagery2)
Nway = size(imagery1);
if length(Nway)>3
    imagery1 = reshape(imagery1,Nway(1),Nway(2),[]);
    imagery2 = reshape(imagery2,Nway(1),Nway(2),[]);
end
psnr = zeros(prod(Nway(3:end)),1);
ssim = psnr;
for i = 1:prod(Nway(3:end))
    A = imagery1(:, :, i);
    B = imagery2(:, :, i);
    psnr(i) = psnr_index(A, B);
    ssim(i) = ssim_index(A, B);
end
psnr = mean(psnr);
ssim = mean(ssim);

rse = norm(imagery1(:)-imagery2(:)) / norm(imagery1(:));

