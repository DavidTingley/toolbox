function [ mean ] = mean_rates( rates , no_cells, xx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 
% exclude interneurons
for i = 1:no_cells
    r = find(rates(i,:,1));
    if length(r) > 200 
       rates(i,:,1) = 0; 
       no_cells = no_cells - 1;
    end
end
% 

leng = length(xx);
leng = floor(leng/5);
        % averaging 5 pos points to 1 (downsizing)     
      
        for i = 1:no_cells
        c = 1;
        d = 5;
        shrink_rate = 1:leng;
            for b = 1: leng
                shrink_rate(b) = sum(rates(i,c:d,1))/5 +1;  %% avoiding 0-1 range
                c = c + 5;
                d = d + 5;
            end
        shrink_rate = (round(shrink_rate)*100)/100;  %% *100/100 decimal adjust
        new_rates(i,:) = shrink_rate;
        end
        %
        
mean = sum(new_rates(:,:,1));

mean = mean/no_cells;


end

