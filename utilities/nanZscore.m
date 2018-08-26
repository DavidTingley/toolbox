function [data] = nanZscore(data)

idx = find(~isnan(data));
bad = find(isnan(data));
data(idx) = zscore(data(idx));
data(bad) = NaN;
