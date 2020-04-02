function [Sigma]=estimnoise(Image,EchoTime)

xdim=size(Image,1);
ydim=size(Image,2);
tdim=size(Image,3);

deb = 1; 
fin = 10;
xi        = zeros(1,1,fin-deb+1);
xi(1,1,:) = EchoTime(tdim-fin+1:tdim-deb+1);
un        = ones(xdim,ydim,fin-deb+1);
xi        = repmat(xi,xdim,ydim,1);
yi        = Image(:,:,tdim-fin+1:tdim-deb+1);

%%
N = size(xi,3);
Sx = sum(xi,3);
Sx2 = sum(xi.^2,3);
Sy = sum(yi,3);
Sxy = sum(xi.*yi,3);
den = (N*Sx2 - Sx.^2);
a = (N*Sxy - Sx.*Sy)./den;
b = (Sx2.*Sy - Sx.*Sxy)./den;

%%
r2 = (yi-repmat(a,1,1,N).*xi - repmat(b,1,1,N).*un).^2;
r1 = (yi-repmat(a,1,1,N).*xi - repmat(b,1,1,N).*un);
m2 = mean(r2,3)-mean(r1,3).^2;

Sigma  = sqrt(2*m2/(4-pi));

