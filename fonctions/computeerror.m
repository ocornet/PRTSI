function [error,bias,sigma,biasp,sigmap,errorp]=computeerror(tetap,tetatest)

diff=(((tetap-tetatest).^2)./tetap.^2);
errorp=sqrt(sum(diff,2)/size(tetap,2))*100;
error=sqrt(sum(diff(:))/(size(tetap,1)*size(tetap,2)))*100;
biasp=mean(tetatest-tetap,2);
bias=mean(abs(biasp./(mean(tetap,2))))*100;
tetahat=mean(tetatest,2);
sigmap=sqrt(mean((tetatest-repmat(tetahat,1,size(tetatest,2))).^2,2));
sigma=mean(abs(sigmap./(mean(tetap,2))))*100;