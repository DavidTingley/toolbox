%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO DO
% fix binning procedure so they overlap instead of border eachother
% find a valid predictability/performance variable


%% create PETHs with varying/overlapping time-bins for all cells
sigs = 'C:\Users\Nitz_Lab\David\Spike_trains\';
% success_markers = 1:20;

load num_cells
cell_num_total = 101
ncells = num_cells;
for rec = 8
    
wire = [];
for cell_num= 1:ncells(rec)

    cell_num_total = cell_num_total+1
if ncells(rec) > 1
    load([sigs 'SIG_' num2str(cell_num_total,'%05d') '.mat'],'spikes_all','tfile_list','success_markers')

    %%%%%%%%%%%%%%%%%%%%  BUILD IN wire seperator
    if cell_num == 1  % only needs to run on the first iteration through recording
        for temp = 1:length(tfile_list)
            wire(temp) = str2num(tfile_list{temp}(5:6));
        end
    end

    %%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(success_markers)
        if length(spikes_all) > success_markers(i,5)*1000
        t = spikes_all(ceil(success_markers(i,1)*1000+.00000000001)+1000:...
            ceil(success_markers(i,1)*1000+.00000000001)+5000);
%        t = spikes_all_control{i}(cell_num,:);
       c = 1;
       for temp = 1:2:4000-2
          tt(c) = sum(t(temp:temp+1));
          c = 1+c; 
       end
       success_trains{i} = tt;
        end
    end

    win = 1;
for i = [1:2:200] % ms window


for t = 1:length(success_trains)

%     s_train_trials(t,:) = smoothts(success_trains{t},'b',i);
        for temp = 1:1999
        if temp>i & temp+i<1999
           s_train_trials(t,temp)=  sum(success_trains{t}(temp-i:temp+i))/(i*2);
        end
        if temp<=i
            s_train_trials(t,temp)=  sum(success_trains{t}(1:temp+i))/(length(1:temp+i));
        end
        if temp+i>=1999
            s_train_trials(t,temp)=  sum(success_trains{t}(temp-i:1999))/(length(temp-i:1999));
        end
        end
        
        s_train_trials_normal(t,:) = success_trains{t};
end
s3{win} = s_train_trials;
clear s_train_trials
win = 1+win;
end
s_all(cell_num,:) = s3;  % cell # X window size X trial # X bins across trial
s_all_norm{cell_num} = s_train_trials_normal;  % cell # X window size X trial # X bins across trial
clear  success* spiral* s3 tt s_train_trials*
end
end

% checks for bad recordings
 if ncells(rec) > 1
bad_rec = 0;
for i = 1:ncells(rec)
    for j =1:ncells(rec)
     if length(s_all{i,1}(:,1)) ~= length(s_all{j,1}(:,1))
        bad_rec = 1; 
     end
    end
end
%% train prediction algorithm
    % HMM, linear regression, GLM
   
        if bad_rec == 0
for cell_num = 1:ncells(rec)
    
    if ~isempty(find(wire(cell_num) ~= wire(:)))
    for win_size = 1:100

            for trial_num = 1:size(s_all{cell_num,win_size},1)

                x = [1:size(s_all{cell_num,win_size},2)];
                y = s_all_norm{cell_num}(trial_num,:)';
                % control using just phase/position of trial

%                 [results_justpos dev_justpos(cell_num,win_size,trial_num)] = glmfit(x',y,'normal');
                [results_justpos] = glmfit(x',y,'normal');
                yhat_justpos(trial_num,:) = glmval(results_justpos,x','identity');         
                
                % both peers and phase/position
                x = [x;];
                for i = 1:ncells(rec)
                    if i ~= cell_num && wire(i) ~= wire(cell_num)
                    x = [x; s_all{i,win_size}(trial_num,:)];
                    end
                end

%                 [results dev(cell_num,win_size,trial_num)] = glmfit(x',y,'normal');
                 [results] = glmfit(x',y,'normal');
                yhat(trial_num,:) = glmval(results,x','identity');
                
                % just peers
                x = [];
                for i = 1:ncells(rec)
                    if i ~= cell_num && wire(i) ~= wire(cell_num)
                    x = [x; s_all{i,win_size}(trial_num,:)];
                    end
                end

%                 [results_justpeers dev_justpeers(cell_num,win_size,trial_num)] = glmfit(x',y,'normal');
                [results_justpeers] = glmfit(x',y,'normal');
                yhat_justpeers(trial_num,:) = glmval(results_justpeers,x','identity');
                
                clear x y
                
                params = [mean(s_all_norm{cell_num}(trial_num,:)),std(s_all_norm{cell_num}(trial_num,:))];
                ll2(cell_num,win_size,trial_num) = normlike(params,s_all_norm{cell_num}(trial_num,:));
                bothll2(cell_num,win_size,trial_num) = normlike(params,mean(yhat(trial_num,:),1));
                peersll2(cell_num,win_size,trial_num) = normlike(params,mean(yhat_justpeers(trial_num,:)));
                posll2(cell_num,win_size,trial_num) = normlike(params,mean(yhat_justpos(trial_num,:)));
                
            end
            
            params = [mean(mean(s_all_norm{cell_num})),std(std(s_all_norm{cell_num}))];
            ll(cell_num,win_size) = normlike(params,mean(s_all_norm{cell_num}));
            bothll(cell_num,win_size) = normlike(params,mean(yhat(:,:)));
            peersll(cell_num,win_size) = normlike(params,mean(yhat_justpeers(:,:)));
            posll(cell_num,win_size) = normlike(params,mean(yhat_justpos(:,:)));  
   
        clear yhat yhat_justpos yhat_justpeers
    end
    end
end

%% test prediction algorithm
% log likelihood


eval(['save C:\Users\Nitz_Lab\David\assemblies\recording_' num2str(rec)])
        end
    end
clearvars -except sigs rec ncells  cell_num_total
end




