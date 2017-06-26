function [ p, place_field ] = calc_place_fields2(i, locs, h, rates_in, field_count, THRESHOLD)
% Calculates place fields

  place_field = 1:length(rates_in(i,:,1));
  place_field(:) = 0;
  y = 0;
  exit = 0;
 
temp_rate = 1;
        p = 1;
        while temp_rate > .1  && exit == 0
            temp_rate = rates_in(i,locs(h)-p,1)/rates_in(i,locs(h),1) + .01;
            p = p + 1;
            if p == locs(h)
                temp_rate = .05;
                p = p - 1;
            end
            if rates_in(i,locs(h) - p,1) < THRESHOLD % minimum firing rate of 5
                for j = 1:10  % Number of position points backward to look 
                    if locs(h) + j <= length(rates_in)
                    if rates_in(i,locs(h)+j,1) > THRESHOLD
                    y = 1;
                    end
                    end
                end
                y
                if y == 0;
                    exit = 1;
                end
            end
        end
temp_rate = 1;
        q = 1;
        while temp_rate > .1 && exit == 0
            temp_rate = rates_in(i,locs(h)+q,1)/rates_in(i,locs(h),1) + .01;
 
            if q + locs(h) == length(rates_in(i,:,1)) 
                temp_rate = .05;
            else 
                q = q + 1;
            end
            if rates_in(i,locs(h) + q,1) < THRESHOLD % minimum firing rate of 5
                for j = 1:10  % Number of position points forward to look 
                    if locs(h) - j > 0
                    if rates_in(i,locs(h)-j,1) > THRESHOLD
                    y = 2;
                    end
                    end
                end
                y
                if y == 0;
                    exit = 2;
                end
            end
        end

        
        if q ~= 1 || p ~= 1 
            if length(locs(h)-p:locs(h)+q) > 8
            place_field(locs(h)-p:locs(h)+q)= 10 + field_count; % rates_in(i,locs(h)-p:locs(h)+q,1);
            end
        else
            place_field(:) = 0;
        end

end

