
function assembly_analysis_result_compiler_resorted_spikes

        d = [];
        dc = [];
        dpeer = [];
        dpos = [];
        list = [];
        index = [];

for batch = 1:length(dir('/home/david/nitzlabdata/assemblies_spike_resort/*'))-2


       for win = 1:300
       e(win) = exist(['/home/david/nitzlabdata/assemblies_spike_resort/' num2str(batch) ...
            '/population_assembly_resort_spikes_winsize_'  num2str(win-1) '.mat']);
       end
       
       if sum(e) == 600
       [dev_new dev_justpeers_new dev_justpos_new dev_control_new cell_list ind indc]...
            = compile_data(batch);
        
%         count = [count;repmat(c,size(cell_list,1),1)];
%         rec_place1 = [rec_place1; sum(num_cells(1:c-1))+cell_list(:,1)];
%         rec_place2 = [rec_place2; sum(num_cells(1:c-1))+cell_list(:,2)];
        list = [list;cell_list];
        d = [d;mean(dev_new,2)];
        dc = [dc;mean(dev_control_new,2)];
        dpeer = [dpeer;mean(dev_justpeers_new,2)];
        dpos = [dpos;mean(dev_justpos_new,2)];
        index = [index;ind'];
%         values = [values;val'];
%         indexc = [indexc;indc'];
%         valuesc = [valuesc;valc'];
        
        
        
        e = zeros(200,1);
 
       end 
end 

   
save assembly_population_results_resorted_spikes.mat
 
end

 
function   [dev_new dev_justpeers_new dev_justpos_new dev_control_new list ind indc] = compile_data(batch) 
  
    for win = 1:300

      load(['/home/david/nitzlabdata/assemblies_spike_resort/' num2str(batch) ...
            '/population_assembly_resort_spikes_winsize_'  num2str(win-1) '.mat']...
            ,'dev*','list')
  
% 
%       ll_new(:,win) = ll;
%       ll2_new(:,:,win) = ll2;
% 
%       posll_new(:,win) = posll;
%       posll2_new(:,:,win) = posll2;
% 
%       bothll_new(:,win) = bothll; 
%       bothll2_new(:,:,win) = bothll2;  
% 
%       peersll_new(:,win) = peersll;
%       peersll2_new(:,:,win) = peersll2;

      dev_new(:,:,win) = squeeze(cell2mat(dev'));
      dev_justpeers_new(:,:,win) = squeeze(cell2mat(dev_justpeers'));
      dev_justpos_new(:,:,win) = squeeze(cell2mat(dev_justpos'));
      dev_control_new(:,:,win) = squeeze(cell2mat(dev_control'));
      
%       for i = 1:size(yhat,1)
%           for j =1:size(yhat,2)
%               MSE_both(i,j,win) = mean((s_all_norm(i,j,:) - yhat(i,j,:)).^2);
%               MSE_peers(i,j,win) = mean((s_all_norm(i,j,:) - yhat_justpeers(i,j,:)).^2);
%               MSE_pos(i,j,win) = mean((s_all_norm(i,j,:) - yhat_justpos(i,j,:)).^2);
%           end
%       end 
      clear ll ll2 posll posll2 peersll peersll2 bothll bothll2 dev dev_justpeers dev_justpos yhat* s_all
    end 

    
        for i = 1:size(dev_new,1)
          [val(i) ind(i) ] = min(mean(squeeze(dev_new(i,:,:))));
          [valc(i) indc(i) ] = min(mean(squeeze(dev_control_new(i,:,:))));
        end
    
    

end


