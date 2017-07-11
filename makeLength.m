function [interped] = makelength(signal,len)
% figure it out.

    interped = interp1(1:length(signal),signal,1:length(signal)/len:length(signal));
c=1;
while length(interped) < len
    
    interped = interp1(1:length(signal),signal,1:length(signal)/(len+c):length(signal));
    
    c = 1+c;
end
end