function [ImageNoise] = addNoise(Image,SNR)
Image1=squeeze(Image(:,:,1));
if(SNR>50)
sd=mean(Image1(Image1~=0))/SNR;
else
    sd=SNR;
end
taillex=size(Image,1);
tailley=size(Image,2);
echonumber=size(Image,3);
phi=2*pi*rand(taillex,tailley,echonumber);
xn = Image.*cos(phi) + sd*randn(size(Image));
yn= Image.*sin(phi) + sd*randn(size(Image));
ImageNoise=sqrt(xn.^2+yn.^2);
