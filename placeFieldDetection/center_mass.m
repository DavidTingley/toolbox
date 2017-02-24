function  centermass = center_mass(gmax_field,rates)
%Calculates the center of mass of all fields

i = length(gmax_field(:,1,1));
j = length(gmax_field(1,1,:));
ind = zeros(i,j);
centermass = zeros(i,j);
first = zeros(i,j);
last = zeros(i,j);
for q = 1:i
    for p = 1:j
    
        first(q,p) = sum(find(gmax_field(q,:,p),1,'first'));
        last(q,p) = sum(find(gmax_field(q,:,p),1,'last'));
        if first(q,p) ~= 0 && last(q,p) ~= 0
        ind(q,p) = sum(rates(q,first(q,p):last(q,p),1));
        end
    end
end         
for q = 1:i
    for p = 1:j
        Mtotal = ind(q,p);
        COM = ind(q,p)/2;
        
        for h = 1:length(rates)
            
           Mtotal = Mtotal - rates(q,first(q,p)+h,1);
          
           if Mtotal < COM
               if Mtotal > 0
                centermass(q,p) = first(q,p)+h; 
                break;
               end
           end 
        end
    end
end


%         e = 0;
%         z = 0;
% %         if i == 5  % cell number with multiple fields
% %         index = 7;  %column number of field_width
% %         end
%         Mtotal = sum(place_field(:,index));
%         COM = Mtotal/2;
%         while z ~= 1;
%         
%             Mtotal = Mtotal - place_field(locs(index)-p+e,index); 
%             e = e + 1;
%             if COM > Mtotal || Mtotal == 0
%                 z = 1;
%             end             
%             
%         end
%         
%         centermass(i,index) = locs(index)-p+e;
end

