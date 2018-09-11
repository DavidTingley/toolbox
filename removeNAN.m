function [data] = removeNAN(data)

data = data(~isnan(data));