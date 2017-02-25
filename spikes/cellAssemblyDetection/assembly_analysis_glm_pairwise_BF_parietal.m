% This file runs the GLM through all available spike trains
function assembly_analysis_glm_pairwise_BF_parietal
% list = dir('//media/nitzlab/Iomega HDD/David/recording_trains_pairwise/');
% list = list(3:end);
% load('/media/nitzlab/Iomega HDD/David/w_parietal.mat','w_parietal');
% w = w_parietal; clear w_parietal;
% load('/media/nitzlab/Iomega HDD/David/spike_overlap_percentages_within_3ms_parietal.mat','list_percentages');

%% should run once to do all glm analyses
rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
load('/media/nitzlab/Iomega HDD/David/parietal_dataset.mat','file','num_cells');
ncells = num_cells;

for rat = [1:8]
for recording = 1:99

       e = exist(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise_parietal/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(1) '.mat']);
       n=zeros(length(file(1).name),1);
       
       for i = 1:length(ncells)
           n = file(i).name==[rats(rat,:) '_rec' sprintf('%02.0f',recording) '_ratemap_anal_parietal.mat']; 
           if sum(n) == length(file(i).name)
              count = i;
              n=zeros(length(file(i).name),1);
           end
       end
       if e ~= 0
           mkdir(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',recording)])
           
           
         p = load(['/media/nitzlab/Iomega HDD/David/new_ratemaps/parietal/' ...
             rats(rat,:) '_rec' sprintf('%02.0f',recording) '_ratemap_anal_parietal.mat']);
         b = load(['/media/nitzlab/Iomega HDD/David/new_ratemaps/' ...
             rats(rat,:) '_rec' sprintf('%02.0f',recording) '_ratemap_anal.mat']);
         for i = 1:size(p.success_timenorm_rates_smooth,1)
             for j =1:size(b.success_timenorm_rates_smooth,1)
                temp = corrcoef(p.success_timenorm_rates_smooth_means_stds_stes(i,:,1),...
                     b.success_timenorm_rates_smooth_means_stds_stes(j,:,1));
                  ww(i,j) = temp(2);
             end
         end
         w = zeros(size(ww,1)+size(ww,2));
         w(1:size(ww,1),size(ww,1)+1:end) = ww;
         w(size(ww,1)+1:end,1:size(ww,1)) = ww';
         
       % splice variables
%        if ncells(count) > 1
%        for win = 0:600
%           ww{win+1} = w(sum(ncells(1:count-1))+1:sum(ncells(1:count)),sum(ncells(1:count-1))+1:sum(ncells(1:count)));
%           list_p{win+1} = list_percentages(sum(ncells(1:count-1))+1:sum(ncells(1:count)),sum(ncells(1:count-1))+1:sum(ncells(1:count)));
%        end
            
       parfor win = 0:600
       if ~exist(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'])
        parGLMRun(win,rat,rats,recording,w);
       end
       end
       clear ww
       end
%        end
end
end
return



function  parGLMRun(win,rat,rats,recording,ww)
%% train prediction algorithm 
%     % HMM, linear regression, GLM  


    cell_pair_num = 0;
    cell_list = [];
    
    par = load(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise_parietal/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'],...
        's_all','s_all_norm');
    
    bf = load(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'],...
        's_all','s_all_norm');
 
    s_all = [par.s_all;bf.s_all];
    s_all_norm = [par.s_all_norm;bf.s_all_norm];
    wire = [ones(size(par.s_all,1),1);ones(size(bf.s_all,1),1).*2];
    clear par bf

%     if length(unique(wire))>1
    for cell_num = 1:length(wire)
    if ~isempty(find(wire(cell_num)==wire(:))) 
    for cnum2 = 1:length(wire)  % should be cell_num or 1?
        if cnum2 ~= cell_num
            if wire(cnum2) ~= wire(cell_num)
            if abs(ww(cell_num,cnum2)) > .6
%             if list_p(cell_num,cnum2) < .2    

%             if list_p(cnum2,cell_num) < .2       
                  cell_pair_num = 1 + cell_pair_num;
           

%     yhat* the actual model fits, and stats* the model stats are TOO LARGE
%     to be saved
%     [yhat{cell_pair_num}, yhat_justpos{cell_pair_num}, yhat_justpeers{cell_pair_num}, yhat_control{cell_pair_num}...
%     stats{cell_pair_num},stats_justpeers{cell_pair_num},stats_justpos{cell_pair_num},stats_control{cell_pair_num}...
         
%     [yhat{cell_pair_num}, yhat_justpos{cell_pair_num}, yhat_justpeers{cell_pair_num}, yhat_control{cell_pair_num}...
        [ll2{cell_pair_num},bothll2{cell_pair_num},posll2{cell_pair_num},peersll2{cell_pair_num},controll2{cell_pair_num}...
    ,dev{cell_pair_num},dev_justpos{cell_pair_num},dev_justpeers{cell_pair_num},dev_control{cell_pair_num},...
     ] = glm_for_trials_pairwise(s_all,s_all_norm,cell_num,cnum2);
    
            cell_list = [cell_list; cell_num cnum2];
                    
%             params = [mean(mean(s_all_norm(cell_num,:,:))),std(std(s_all_norm(cell_num,:,:)))];
%             ll{cell_pair_num} = normlike(params,mean(s_all_norm(cell_num,:,:)));
%             bothll{cell_pair_num} = normlike(params,mean(yhat{cell_pair_num}));
%             peersll{cell_pair_num} = normlike(params,mean(yhat_justpeers{cell_pair_num}));
%             posll{cell_pair_num} = normlike(params,mean(yhat_justpos{cell_pair_num}));  
   
            end
            end
%             end
%             end
        end
    end
    end
    end

if ~isempty(cell_list)
        disp(['file ' rats(rat,:) '_rec_' sprintf('%02.0f',recording)...
            '_winsize_'  num2str(win) '.mat saved...'])
       
        save(['/media/nitzlab/Iomega HDD/David/assemblies_cross/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'],...
           '*control*','dev*','cnum2','wire','*list')%,'*ll*','params','*yhat*') % *stats*
end
%     end

return