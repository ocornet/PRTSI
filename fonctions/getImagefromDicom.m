function [Image]=getImagefromDicom(path) 
fclose all;
list=dir([path '*.IMA']);
for k=1:(length(list))
Imageze2(:,:,k) = dicomread([path list(k).name]);
end
Image=double(Imageze2);
