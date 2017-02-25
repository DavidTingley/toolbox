%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO DO
% fix binning procedure so they overlap instead of border eachother
% find a valid predictability/performance variable

%% create PETHs with varying/overlapping time-bins for all cells
sigs = '/home/nitzlab/Dropbox/Spike_trains/';

load /home/nitzlab/Dropbox/num_cells
ncells = num_cells;
% if rec_start > 1 
%     cell_num_total = sum(ncells(1:rec_start-1));
% elseif rec_start == 1
%     cell_num_total = 0;
% end

% 
% for rec = rec_start
% 
% wire = [];
% for cell_num= 1:ncells(rec)
% 
%     cell_num_total = cell_num_total+1
% if ncells(rec) > 1
% 
%     % bad_rec = 0;
% % for i = 1:ncells(rec)
% %     for j =1:ncells(rec)
% %      if length(s_all{i}(1,:,1)) ~= length(s_all{j}(1,:,1))
% %         bad_rec = 1; 
% %      end
% %     end
% % end
% 
%     load([sigs 'SIG_' num2str(cell_num_total,'%05d') '.mat'],'spikes_all',...
%         'tfile_list','success_markers','recording','ratname')
%     %%%%%%%%%%%%%%%%%%%%  BUILD IN wire seperator
% 
%     if cell_num == 1  % only needs to run on the first iteration through recording
%         wire = zeros(length(tfile_list),1);
%         parfor t= 1:length(tfile_list)
%             wire(t) = str2double(tfile_list{t}(5:6));
%         end
%     end 
% 
%     %%%%%%%%%%%%%%%%%%%%
%     success_trains = spalloc(length(success_markers),1999,.02*4000*length(success_markers));
%     parfor i = 1:length(success_markers)
%          success_trains(i,:) = sparse(test_func2(success_markers,spikes_all,i));
%     end
%     
% tic
% parfor i = 1:200 % ms window
%       s_all= test_func(success_trains,i);
%       s_all_norm = success_trains;
        save(['/home/nitzlab/recording_trains/' ratname '_' recording '_winsize_' i],'-v7.3');
% end
% toc
% clear  success* spiral* tt
% end 
% end

% clearvars -except rec sigs cell_num_total ncells temp
% 
% end



list = dir('/home/nitzlab/recording_trains/');
list = list(3:end);

for rec = 1
  tic  
load(['/home/nitzlab/recording_trains/' list(rec).name],'s_all','s_all_norm','wire')

bad_rec = 0;
for i = 1:ncells(rec)
    for j =1:ncells(rec)
     if length(s_all{i}(1,:,1)) ~= length(s_all{j}(1,:,1))
        bad_rec = 1; 
     end
    end
end

%% train prediction algorithm
%     % HMM, linear regression, GLM
    if ncells(rec) > 1
        if bad_rec == 0
for cell_num = 1:ncells(rec)
    if length(unique(wire))>1
        %% can we slice s_all and use a parfor loop?
    for win_size = 1:200

    [yhat, yhat_justpos, yhat_justpeers, ll2,bothll2,posll2,peersll2,dev(cell_num,win_size,:)...
        ,dev_justpos(cell_num,win_size,:),dev_justpeers(cell_num,win_size,:)...
        ] = glm_for_trials(s_all,s_all_norm,cell_num,win_size,wire);
            
            params = [mean(mean(s_all_norm{cell_num})),std(std(s_all_norm{cell_num}))];
            ll(cell_num,win_size) = normlike(params,mean(s_all_norm{cell_num}));
            bothll(cell_num,win_size) = normlike(params,mean(yhat(:,:)));
            peersll(cell_num,win_size) = normlike(params,mean(yhat_justpeers(:,:)));
            posll(cell_num,win_size) = normlike(params,mean(yhat_justpos(:,:)));  
   
        
    end
    end
disp(cell_num)
end

%% test prediction algorithm
% log likelihood

clear s_all s_all_norm
time = toc
eval(['save /home/nitzlab/assemblies/' list(rec).name])
        end
    end
clearvars -except sigs rec ncells  cell_num_total

end




