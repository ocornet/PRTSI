function mask = mymask(im,ths)
dimx = size(im,1);
dimy = size(im,2);

mask = zeros(dimx,dimy);
for k=1:dimx
    for l=1:dimy
        mask(k,l)=((l-1)*dimx)+k;
    end
end
mask = mask(sum(im,3)>ths*max(max(sum(im,3))));