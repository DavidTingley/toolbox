
function assembly_analysis_result_compiler

rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];


for rat =   1:8
   parfor rec = 1:100 % can be a parfor loop if glm_test isnt running 
       for win = 1:200
       e(win) = exist(['/media/nitzlab/Iomega HDD/David/recording_trains_smart/' rats(rat,:) '_rec_'  sprintf('%02.0f',rec)...
           '_winsize_'  num2str(win) '.mat']);
       end
       
       if sum(e) == 400
        compile_data(rats,rat,rec)
        e = zeros(200,1);
 
       end 
   end 
end  
end 
 
  
 
function  compile_data(rats,rat,rec) 
  
    for win = 1:200

      load(['/media/nitzlab/Iomega HDD/David/recording_trains_smart/' rats(rat,:) '_rec_'  sprintf('%02.0f',rec)...
       '_winsize_' num2str(win) '.mat'],'both*','peers*','pos*','ll*','dev*','yhat*','s_all_norm','wire')
    if length(unique(wire))>1

      ll_new(:,win) = ll;
      ll2_new(:,:,win) = ll2;

      posll_new(:,win) = posll;
      posll2_new(:,:,win) = posll2;

      bothll_new(:,win) = bothll; 
      bothll2_new(:,:,win) = bothll2;  

      peersll_new(:,win) = peersll;
      peersll2_new(:,:,win) = peersll2;

      dev_new(:,:,win) = dev;
      dev_justpeers_new(:,:,win) = dev_justpeers;
      dev_justpos_new(:,:,win) = dev_justpos;

      for i = 1:size(yhat,1)
          for j =1:size(yhat,2)
              MSE_both(i,j,win) = mean((s_all_norm(i,j,:) - yhat(i,j,:)).^2);
              MSE_peers(i,j,win) = mean((s_all_norm(i,j,:) - yhat_justpeers(i,j,:)).^2);
              MSE_pos(i,j,win) = mean((s_all_norm(i,j,:) - yhat_justpos(i,j,:)).^2);
          end
      end 
      clear ll ll2 posll posll2 peersll peersll2 bothll bothll2 dev dev_justpeers dev_justpos yhat* s_all
    end 
    end

if length(unique(wire))>1
   
      disp([rats(rat,:) '_rec_'  sprintf('%02.0f',rec) ' file done...'])
      save(['/media/nitzlab/Iomega HDD/David/assemblies/' rats(rat,:) '_rec_' sprintf('%02.0f',rec)])        

end
end


