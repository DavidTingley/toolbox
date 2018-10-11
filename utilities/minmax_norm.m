function [m] = minmax_norm(mm)

m = mm - min(mm);
m = m ./ max(m);
return 