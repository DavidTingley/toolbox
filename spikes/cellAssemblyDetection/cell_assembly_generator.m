for ii = 1:20
    
    clearvars -except ii spikes_all*
    NAC_Fundamentals
    
   t = zeros(42,6000);
for i = 1:21
for j = 1:6000
for f = 1:length(Spks)
if Spks(f,1) == j & Spks(f,2) == i
t(i,j) = 1;
end
end
end
end 

% r = 40 + randi(10,1);
% rr = 40 + randi(10,1);
spikes_all(ii,:,:) = t; %[t(1,:);[t(1,r:end),zeros(1,r-1)];[t(1,rr:end),zeros(1,rr-1)]];



end


%% control set

% for ii = 1:20
%     
%     clearvars -except ii spikes_all* 
% 
% 
%    for control = 1:21
%     NAC_Fundamentals
%     t = zeros(21,6000);
%     r = randperm(21);
% for i = 1:21
% for j = 1:6000
% for f = 1:length(Spks)
% if Spks(f,1) == j & Spks(f,2) == i
% t(i,j) = 1;
% end
% end
% end
% end
% train(control,:) = t(r(control),:);
%    end
% spikes_all_control{ii} = train;  % should have no cell assembly structure
% 
% end
