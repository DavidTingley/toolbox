function [dist] = fracdist(x,y,frac)

dist = sum((x - y).^frac).^(1/frac);

return