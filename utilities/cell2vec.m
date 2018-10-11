function [vec] = cell2vec(cell)

vec =[];
for i = 1:length(cell)
    vec = [vec;cell{i}(:)];
end