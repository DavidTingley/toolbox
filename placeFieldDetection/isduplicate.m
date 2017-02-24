function [cell_field] = isduplicate(cell_field, no_cells, xx, ...
    max_loc, no_maxs, index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
 
for i = 1:no_cells
        dup = 1:no_maxs(i);
        dup(:) = 0;
 
        for j = 1:no_maxs(i)
            ind = find(cell_field(i,:));
            cell_field(i,ind(:)) = 10;            

            emp = find(cell_field(i,:,j));
        if emp ~= 0
                                 
            for k = 1:length(xx)
                if cell_field(i,k,j) == cell_field(i,k,index(i))
                    if cell_field(i,k,j) ~= 0
                    dup(j) = 1;
                    end
                end
            end
        end
        end 
        for g = 1:no_maxs(i)
             if dup(g) == 1 && g ~= index(i)
                 cell_field(i,:,g) = 0;
             end
        end       
        
       %% now I need to deal with seperate fields that dont overlap 
        no_overlap = 1:no_maxs;
        no_overlap(:) = 0;
        
       for m = 1:no_maxs(i)
           if m ~= index(i)
               for n = 1:no_maxs(i)
                   if m ~= n
                       for k = 1:length(xx)
                           if cell_field(i,k,m) == cell_field(i,k,n)
                               if cell_field(i,k,m) ~= 0
                                   dup(m) = 2;
                               end
                           end
                       end
                   end
               end     
           end 
        end
               
           %% delete zeros and none-global max duplicates
           for w = 1:no_maxs(i)
            if dup(w) == 2
                  no_overlap(w) = length(find(cell_field(i,:,w))); 
            end
           end
           
           [C I] = max(no_overlap);
           for m = 1:no_maxs(i)
               if dup(m) == 2
                  if I ~= m
                     cell_field(i,:,m) = 0;
 
                  end
               end
           end
end
end
