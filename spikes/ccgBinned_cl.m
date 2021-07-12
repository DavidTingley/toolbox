function [ccg pvals] = ccgBinned_cl(sig1,sig2,nBins)

sig1 = squeeze(sig1)';
sig2 = squeeze(sig2)';

for i=-nBins:nBins
    [ccg(i+nBins+1) pvals(i+nBins+1)] = circ_corrcl(circshift(sig1,i),sig2);
end

