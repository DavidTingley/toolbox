function [yhat, yhat_justpos, yhat_justpeers, ll2,bothll2,posll2,peersll2...
    dev,dev_justpos,dev_justpeers]= glm_for_trials(s_all,s_all_norm,cell_num,win_size,wire)

yhat_justpos = zeros(size(s_all{1}{1}));
yhat_justpeers = zeros(size(s_all{1}{1}));
yhat = zeros(size(s_all{1}{1}));

for trial_num = 1:size(s_all{cell_num}{1},1)

                x = [1:size(s_all{cell_num}{win_size},2)];
                y = s_all_norm{cell_num}(trial_num,:)';
                % control using just phase/position of trial

                [results_justpos dev_justpos(trial_num)] = glmfit(x',y,'normal');
%                 [results_justpos dev_justpos] = glmfit(x',y,'normal');
                yhat_justpos(trial_num,:) = glmval(results_justpos,x','identity');         
                
                % both peers and phase/position
                x = [x;];
                for i = 1:length(s_all)
                    if i ~= cell_num && wire(i) ~= wire(cell_num)
                    x = [x; s_all{i}{win_size}(trial_num,:)];
                    end
                end
                
                [results dev(trial_num)] = glmfit(x',y,'normal');
%                  [results dev] = glmfit(x',y,'normal');
                yhat(trial_num,:) = glmval(results,x','identity');
                
                % just peers
                x = [];
                for i = 1:length(s_all)
                    if i ~= cell_num && wire(i) ~= wire(cell_num)
                    x = [x; s_all{i}{win_size}(trial_num,:)];
                    end
                end
                if ~isempty(x)
                [results_justpeers dev_justpeers(trial_num)] = glmfit(x',y,'normal');
%                 [results_justpeers dev_justpeers] = glmfit(x',y,'normal');
                yhat_justpeers(trial_num,:) = glmval(results_justpeers,x','identity');
                end


                
                params = [mean(s_all_norm{cell_num}(trial_num,:)),std(s_all_norm{cell_num}(trial_num,:))];
                ll2(cell_num,win_size,trial_num) = full(normlike(params,s_all_norm{cell_num}(trial_num,:)));
                bothll2(cell_num,win_size,trial_num) = full(normlike(params,(yhat(trial_num,:))));
                peersll2(cell_num,win_size,trial_num) = full(normlike(params,(yhat_justpeers(trial_num,:))));
                posll2(cell_num,win_size,trial_num) = full(normlike(params,(yhat_justpos(trial_num,:))));
                
end
            
end