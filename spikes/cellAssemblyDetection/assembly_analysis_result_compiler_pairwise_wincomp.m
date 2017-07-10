% This funciton compiles GLM results into a more readable format
function assembly_analysis_result_compiler

rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];


for rat = 1:8
   parfor rec = 1:100 % can be a parfor loop if glm_test isnt running 
       for win = 1:601
       e(win) = exist(['F:\David\assemblies\' rats(rat,:) '\' sprintf('%02.0f',rec) '\'...
           rats(rat,:) '_rec_'  sprintf('%02.0f',rec) '_winsize_'  num2str(win-1) '.mat']);
       end
       if sum(e) == 1200
        compile_data(rats,rat,rec)
        e = zeros(600,1);
       end
   end
end
end
  
 
 
function  compile_data(rats,rat,rec)
  
    for win = 1:601

      load(['F:\David\assemblies\' rats(rat,:) '\' sprintf('%02.0f',rec) '\'...
           rats(rat,:) '_rec_'  sprintf('%02.0f',rec) '_winsize_'  num2str(win-1) '.mat'],'both*',...
           'peers*','pos*','ll*','control*','dev*','yhat*','s_all_norm','wire','*list*','stats*')
    if length(unique(wire))>1
      ll_new(:,win) = cell2mat(ll);
      ll2_new(:,:,win) = cell2mat(ll2);

      posll_new(:,win) = cell2mat(posll);
      posll2_new(:,:,win) = cell2mat(posll2);

      bothll_new(:,win) = cell2mat(bothll);
      bothll2_new(:,:,win) = cell2mat(bothll2);

      peersll_new(:,win) = cell2mat(peersll);
      peersll2_new(:,:,win) = cell2mat(peersll2);

      dev_new(:,:,win) = dev;
      dev_justpeers_new(:,:,win) = dev_justpeers; 
      dev_justpos_new(:,:,win) = dev_justpos;

%       for i = 1:size(yhat,2)
%           for j =1:size(yhat{1},1)
%               MSE_both(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat{i}(j,:)).^2;
%               MSE_peers(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat_justpeers{i}(j,:)).^2;
%               MSE_pos(i,j,win) = mean(squeeze(s_all_norm(i,j,:))' - yhat_justpos{i}(j,:)).^2;
%           end
%       end 
      clear ll ll2 posll posll2 peersll peersll2 bothll bothll2 dev dev_justpeers dev_justpos yhat* s_all*
    end 
    end

if length(unique(wire))>1
      disp([rats(rat,:) '_rec_'  sprintf('%02.0f',rec) ' file done...'])
      save(['F:\David\results\'...
           rats(rat,:) '_rec_'  sprintf('%02.0f',rec)  '.mat'])        
end
end


