function [interped] = makelength(signal,len)
% figure it out.
interped = interp1(1:length(signal),signal,1:length(signal)/len:length(signal));

end