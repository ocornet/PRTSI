function I = unmask(masked,dimx,dimy,nbComponents, mask)

I = zeros(dimx, dimy, nbComponents);
for i=1:length(mask)
    [xi,yi] = myunravel(mask(i),dimx);
    I(xi,yi,:) = masked(i,:);
end