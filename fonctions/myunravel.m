function [x,y] = myunravel(z,dimx)
    y = floor((z-1)./dimx)+1;
    x = mod((z-1), dimx)+1;
    
    