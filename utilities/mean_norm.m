function [m] = mean_norm(mm)

for i = 1:size(mm,1)
    m(i,:) = mm(i,:) ./ (nanmean(mm(i,:))+.00000001);
end

return 