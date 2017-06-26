function [cell_field, index, center_m] = ...
    place_field_info(rates, in_template)
% Written by David Tingley for the Nitz Laboratory at UCSD
% 3/17/11    
% This program classifies a place 'field'
% by its length, firing rate, and consistency
    
%    USEFULL VARIABLES
% center_m(no_cells) = the center of mass for the global max's 
%                      field for each cell
% no_maxs(i) = number of maximums for each cell
% cell_field(i, location/length of track, maximum #) =...
%                                   array containing the
%                                   location and length of 
%                                   each identified field
% index(i) = the locations of all of the maxima for each cell
% gmax_field = "field" associated with the global max for the cell
% field_width(no_cells,2) = column 1 

% INPUT
%
% OUTPUT
%
% TO DO
%   - Include parts of field after dip below 10%
%   - NEW FUNCITON = mean firing rate across all cells
%                  = discritize first, and seperate interneurons
%                  = 
%   - Change centermass function to calc multiple fields
%                  = same with field_widths
 
no_cells=size(rates,1);
THRESHOLD = 5; % Minimum firing rate in hertz

   no_maxs = zeros(no_cells,1);
   max_loc = zeros(no_cells,1);
   index = zeros(no_cells,1);
   
for i = 1:no_cells
       field_count = 0;
        % smoothing
        %rates(i,:,1) = smoothts(rates(i,:,1),'g', 10, 10);
        %

     [pks, locs] = findpeaks(rates(i,:,1));
        
    L = length(locs);
    
    [C ind] = max(pks);
    no_maxs(i) = L;
    if L == 0
        index(i) = 0;
        max_loc(i) = 0;
    else
        max_loc(i) = locs(ind);
        index(i) =ind;
    end
        
        


    cell_field2 = zeros(no_cells,length(rates),max(no_maxs));
           

    place_field = zeros(length(rates),L);

    if L == 0
        cell_field(i,:,:) = cell_field2(i,:,:);
        center_m(i) = 0;
    end
    for h = 1:L
    
        [pks, locs] = findpeaks(rates(i,:,1));
        
        L = length(locs);
    
        [p, field] = calc_place_fields2(i, locs, ...
                h, rates, field_count, THRESHOLD);      
            
        field_count = field_count +1;
            
        place_field(:,h)= field(:);
        cell_field(i,:,h) = place_field(:,h);

    end
       
end


% checks for and removes duplicate fields
   [gmax_field] = isduplicate(cell_field, ...
      no_cells, in_template(2:end-1,3), max_loc, no_maxs, index);
  
   [mean_rate] = mean_rates(rates, no_cells, in_template(2:end-1,3));
 
   center_m = center_mass(gmax_field,rates);
  
   field_width = zeros(no_cells,max(no_maxs));

x = 1:length(in_template(2:end-1,3));
 for i = 1:no_cells 

     for j = 1:no_maxs(i)         
    field_width(i,j) = length(find(gmax_field(i,:,j)));
     end
     
        figure(1)
        hold on
        bar(rates(i,:,1)); 
   
        
    for t = 1:no_maxs(i)
        if index(i) ~= 0
        plot(x,gmax_field(i,x,t),'.m');  % global maximum
        end
%        plot(x,cell_field(i,x,t),'.g');     % local maxima
        
    end
    pause; clf;
 end

end

