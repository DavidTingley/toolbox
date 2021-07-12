% This funciton compiles GLM results into a more readable format
function assembly_analysis_result_compiler_pairwise

rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
s = load('/media/nitzlab/data/Dropbox/Dropbox/parietal_dataset.mat','file','num_cells');
file = s.file; num_cells = s.num_cells;
disp('starting to compile data')
warning off

        d = [];
        dc = [];
        dpeer = [];
        dpos = [];
        index = [];
        values = [];
        indexc = [];
        valuesc = [];
        list = [];
        count = [];
        rec_place1 = [];
        rec_place2 = [];  
        trials_dev_all = [];
        trials_dev_control_all = [];
        
for rat = 1:8
   for rec = 1:99 
       
        for i = 1:length(num_cells)  % can remove this once count is saved by assembly_analysis_glm_pairwise
           n = file(i).name==[rats(rat,:) '_rec' sprintf('%02.0f',rec) '_ratemap_anal_parietal.mat']; 
           if sum(n) == length(file(i).name)
              c = i;
              n=zeros(length(file(i).name),1);
           end
        end
       
        
       if exist(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',rec) ])
        e = dir(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',rec) '/'...
            rats(rat,:) '_rec_'  sprintf('%02.0f',rec) '_winsize_*']);
       if length(e) == 601
           
        [trials_dev,trials_dev_control,dev_new dev_justpeers_new dev_justpos_new dev_control_new ind val indc valc cell_list]...
            = compile_data(rats,rat,rec);
        
        trials_dev_all = [trials_dev_all; trials_dev];
        trials_dev_control_all = [trials_dev_control_all; trials_dev_control];
        
        count = [count;repmat(c,size(cell_list,1),1)];
        rec_place1 = [rec_place1; sum(num_cells(1:c-1))+cell_list(:,1)];
        rec_place2 = [rec_place2; sum(num_cells(1:c-1))+cell_list(:,2)];
        list = [list;cell_list];
        d = [d;mean(dev_new,3)];
        dc = [dc;mean(dev_control_new,3)];
        dpeer = [dpeer;mean(dev_justpeers_new,3)];
        dpos = [dpos;mean(dev_justpos_new,3)];
        index = [index;ind'];
        values = [values;val'];
        indexc = [indexc;indc'];
        valuesc = [valuesc;valc'];
       end
       end
   end
end

save assembly_population_results_cross_region
end
  
 
 
function   [trials_dev,trials_dev_control,dev_new dev_justpeers_new dev_justpos_new dev_control_new ind val indc valc cell_list]...
    = compile_data(rats,rat,rec)
  
    for win = 1:601

      load(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',rec) '/'...
           rats(rat,:) '_rec_'  sprintf('%02.0f',rec) '_winsize_'  num2str(win-1) '.mat']) %,'both*',...
           %'peers*','pos*','ll*','control*','dev*','wire','*list*','count') %,'yhat*','s_all_norm') % *stats*

    if length(unique(wire))>1
%       ll_new(:,win) = cell2mat(ll);
%       ll2_new(:,:,win) = cell2mat(ll2);
% 
%       posll_new(:,win) = cell2mat(posll);
%       posll2_new(:,:,win) = cell2mat(posll2);
% 
%       bothll_new(:,win) = cell2mat(bothll);
%       bothll2_new(:,:,win) = cell2mat(bothll2);
% 
%       peersll_new(:,win) = cell2mat(peersll);
%       peersll2_new(:,:,win) = cell2mat(peersll2);

 for i = 1:size(dev,2)
      trials_dev{i,win} = dev{i};
      trials_dev_control{i,win} = dev_control{i};
            
      dev_new(i,win,:) = mean(dev{i});
      dev_justpeers_new(i,win,:) = mean(dev_justpeers{i}); 
      dev_justpos_new(i,win,:) = mean(dev_justpos{i});
      dev_control_new(i,win,:) = mean(dev_control{i});
 end
%       for i = 1:size(yhat,2)
%           for j =1:size(yhat{1},1)
%               MSE_both(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat{i}(j,:)).^2;
%               MSE_peers(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat_justpeers{i}(j,:)).^2;
%               MSE_pos(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat_justpos{i}(j,:)).^2;
%           end
%       end 
      clear ll ll2 posll posll2 peersll peersll2 bothll bothll2 dev dev_control dev_justpeers dev_justpos yhat* s_all*
    end 
    end
    
    for i = 1:size(dev_new,1)
          [val(i) ind(i) ] = min(mean(dev_new(i,:,:),3));
          [valc(i) indc(i) ] = min(mean(dev_control_new(i,:,:),3));
    end
    

if length(unique(wire))>1
      disp([rats(rat,:) '_rec_'  sprintf('%02.0f',rec) ' file done...'])
      save(['/media/nitzlab/data/Dropbox/Dropbox/results_cross/'...
           rats(rat,:) '_rec_'  sprintf('%02.0f',rec)  '.mat'])        
end
end


