% [yhat, yhat_justpos, yhat_justpeers, yhat_control,...
%     stats,stats_justpeers,stats_justpos,stats_control,...    
    
% function [yhat, yhat_justpos, yhat_justpeers, yhat_control,...
function [ ll2,bothll2,posll2,peersll2, controll2,dev,dev_justpos,dev_justpeers,dev_control...
    ]= glm_for_trials_pairwise_cross(s_all,s_all_norm,cell_num,cnum2)


%% find trials where the cell doesnt fire
marks = zeros(size(s_all,2),1);
for trial_num = 1:size(s_all,2)
   if length(find(s_all_norm(cell_num,trial_num,:)==1)) < 1
       marks(trial_num) = 1;
   end
   if length(find(s_all_norm(cnum2,trial_num,:)==1)) < 1
       marks(trial_num) = 1;
   end
end
s_all= s_all(:,marks==0,:);
s_all_norm = s_all_norm(:,marks==0,:);
%% alters recording to only analyze trials with spikes

% yhat_justpos = zeros(size(squeeze(s_all(1,:,:))));
% yhat_justpeers = zeros(size(squeeze(s_all(1,:,:))));
% yhat = zeros(size(squeeze(s_all(1,:,:))));
r = randperm(size(s_all,2));
rr = randperm(size(s_all,2)); 

for trial_num = 1:size(s_all,2)

                xpos = [repmat(1:size(s_all,3),length(cnum2),1)];
                y = squeeze(s_all_norm(cell_num,trial_num,:))';
                % control using just phase/position of trial

                [results_justpos dev_justpos(trial_num) stats_justpos(trial_num)]...
                    = glmfit(xpos',y,'normal');
%                   yhat_justpos(trial_num,:) = glmval(results_justpos,xpos','identity');         
                      yhat_justpos = glmval(results_justpos,xpos','identity');  
                
                %% both peers and phase/position
                xboth = [repmat(1:size(s_all,3),length(cnum2),1)];
                xboth = [xboth; squeeze(s_all(cnum2,trial_num,:))];

                
                [results dev(trial_num) stats(trial_num)]...
                    = glmfit(xboth',y,'normal');
%                 yhat(trial_num,:) = glmval(results,xboth','identity');
                yhat = glmval(results,xboth','identity');

                %% just peers
                xpeers = [repmat(squeeze(s_all(cnum2,trial_num,:)),length(cnum2),1)];

                [results_justpeers dev_justpeers(trial_num) stats_justpeers(trial_num)]...
                    = glmfit(xpeers',y,'normal');
%                 yhat_justpeers(trial_num,:) = glmval(results_justpeers,xpeers','identity');
                yhat_justpeers = glmval(results_justpeers,xpeers','identity');
                %% TODO add control (spike trains from another rat/recording)
                

                xcontrol = [repmat(squeeze(s_all(cnum2,r(trial_num),:)),length(cnum2),1)];
                ycontrol = squeeze(s_all_norm(cell_num,rr(trial_num),:))';
                
                [results_control dev_control(trial_num) stats_control(trial_num)]...
                    = glmfit(xcontrol',ycontrol,'normal');
%                 yhat_control(trial_num,:) = glmval(results_control,xcontrol','identity');
                 yhat_control = glmval(results_control,xcontrol','identity');
                
                

                params = [mean(s_all_norm(cell_num,trial_num,:)),std(s_all_norm(cell_num,trial_num,:))];
                
%                 ll2(trial_num) = full(normlike(params,s_all_norm(cell_num,trial_num,:)));
%                 bothll2(trial_num) = full(normlike(params,(yhat(trial_num,:))));
%                 peersll2(trial_num) = full(normlike(params,(yhat_justpeers(trial_num,:))));
%                 posll2(trial_num) = full(normlike(params,(yhat_justpos(trial_num,:))));
%                 controll2(trial_num) = full(normlike(params,(yhat_control(trial_num,:))));
                ll2(trial_num) = full(normlike(params,s_all_norm(cell_num,trial_num,:)));
                bothll2(trial_num) = full(normlike(params,(yhat)));
                peersll2(trial_num) = full(normlike(params,(yhat_justpeers)));
                posll2(trial_num) = full(normlike(params,(yhat_justpos)));
                controll2(trial_num) = full(normlike(params,(yhat_control)));
              
end

            
end