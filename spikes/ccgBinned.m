function [ccg pvals] = ccgBinned(sig1,sig2,nBins)

sig1 = squeeze(sig1)';
sig2 = squeeze(sig2)';

for i=-nBins:nBins
    [ccg(i+nBins+1) pvals(i+nBins+1)] = corr(circshift(sig1,i),sig2,'rows','complete');
end

