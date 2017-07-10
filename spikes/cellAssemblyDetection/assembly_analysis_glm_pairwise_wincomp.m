% This file runs the GLM through all available spike trains
function assembly_analysis_glm_pairwise
% list = dir('\\NITZLAB-SX2801\nitzlab-iomegaHDDavid/recording_trains_pairwise/');
% list = list(3:end);
load('C:\Users\Nitz Lab 2\Dropbox_tingley\Dropbox/w.mat','w');
load('C:\Users\Nitz Lab 2\Dropbox_tingley\Dropbox/spike_overlap_percentages.mat','list_percentages');

%% should run once to do all glm analyses
rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
load('C:\Users\Nitz Lab 2\Dropbox_tingley\Dropbox/ratemaps_10012012_bf.mat','file','num_cells');
ncells = num_cells;

for rat = [6]
for recording = 12%1:99

       e = exist(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/recording_trains_pairwise/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(1) '.mat']);
       n=zeros(length(file(1).name),1);
       
       for i = 1:length(ncells)
           n = file(i).name==[rats(rat,:) '_rec' sprintf('%02.0f',recording) '_ratemap_anal.mat']; 
           if sum(n) == length(file(i).name)
              count = i;
              n=zeros(length(file(i).name),1);
           end
       end
       if e ~= 0
           mkdir(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/assemblies/' rats(rat,:) '/' sprintf('%02.0f',recording)])
           
       % splice variables
       if ncells(count) > 1
       for win = 0:600
          ww{win+1} = w(sum(ncells(1:count-1))+1:sum(ncells(1:count)),sum(ncells(1:count-1))+1:sum(ncells(1:count)));
          list_p{win+1} = list_percentages(sum(ncells(1:count-1))+1:sum(ncells(1:count)),sum(ncells(1:count-1))+1:sum(ncells(1:count)));
       end
           
       parfor win = 400:600
%            if ~exist(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/assemblies/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
%            rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'])
t = dir(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/recording_trains_pairwise/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat']);
       if str2num(t.date(1:2))<15
        parGLMRun(win,rat,rats,recording,list_p{win+1},ww{win+1});
%            end
       end
       end
       end
       end
end
end
return

function  parGLMRun(win,rat,rats,recording,list_p,ww)
%% train prediction algorithm 
%     % HMM, linear regression, GLM  


    cell_pair_num = 0;
    cell_list = [];
    
    load(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/recording_trains_pairwise/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'],...
        's_all','s_all_norm','wire')
 

    if length(unique(wire))>1
    for cell_num = 1:length(wire)
    if ~isempty(find(wire(cell_num)==wire(:))) 
    for cnum2 = 1:length(wire)  % should be cell_num or 1?
        if cnum2 ~= cell_num
            if wire(cnum2) ~= wire(cell_num)
            if abs(ww(cell_num,cnum2)) > .6
            if list_p(cell_num,cnum2) < .2    
            if list_p(cnum2,cell_num) < .2       
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
            end
            end
        end
    end
    end
    end

if ~isempty(cell_list)
        disp(['file ' rats(rat,:) '_rec_' sprintf('%02.0f',recording)...
            '_winsize_'  num2str(win) '.mat saved...'])
       
        save(['\\NITZLAB-SX2801\nitzlab-iomegaHD/David/assemblies/' rats(rat,:) '/' sprintf('%02.0f',recording)  '/'...
           rats(rat,:) '_rec_' sprintf('%02.0f',recording) '_winsize_'  num2str(win) '.mat'],...
           '*control*','dev*','cnum2','wire','*list')%,'*ll*','params','*yhat*') % *stats*
end
    end

return