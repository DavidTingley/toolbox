
% this code takes a folder where individual spike trains are saved and runs
% a pairwise cross correlation on all cell pairs within the folder

list = dir('C:\Users\Nitz_Lab\David\Spike_trains_parietal\*');
load C:\Users\Nitz_Lab\David\ratemaps_09262012_parietal_with_light_offsets
trainpath = 'C:\Users\Nitz_Lab\David\Spike_trains_parietal\';
clearvars -except list trainpath num_cells xcf
tic
list = list(3:length(list));
   for j = 1:sum(num_cells)
      for h = j:sum(num_cells)
          for i = 1:length(num_cells)
              if sum(num_cells(1:i-1)) < j && sum(num_cells(1:i-1)) < h && ...
                      sum(num_cells(1:i)) >= j && sum(num_cells(1:i)) >= h
                 if j~=h
           spikes1 = load([trainpath list(j).name]);
           spikes2 = load([trainpath list(h).name]);
           [xcf_all{j,h},lags,bounds] = crosscorr(spikes1.spikes_all,spikes2.spikes_all,500);
           [xcf_nontrial{j,h},lags,bounds] = crosscorr(spikes1.spikes_nontrial,spikes2.spikes_nontrial,500);
           [xcf_success{j,h},lags,bounds] = crosscorr(spikes1.spikes_success,spikes2.spikes_success,500);
%            [xcf_failure{j,h},lags,bounds] = crosscorr(spikes1.spikes_failure,spikes2.spikes_failure,500);
%            [xcf_nose{j,h},lags,bounds] = crosscorr(spikes1.spikes_nose,spikes2.spikes_nose,500);
           
            j 
            h
                 end
           clearvars -except i j h list trainpath num_cells xcf*           
              end
          end
      end
      toc
   end
   save(['C:\Users\Nitz_Lab\David\crosscor\' 'xcf_parietal.mat']) 
        count = 1;
      for i = 1:length(xcf_all)-1
      for f = 1:length(xcf_all)
          if ~isempty(xcf_success{i,f})
          all_xcorrs(count,:) = xcf_all{i,f};
          success_xcorrs(count,:) = xcf_success{i,f};
%             failure_xcorrs(count,:) = xcf_failure{i,f};
          nontrial_xcorrs(count,:) = xcf_nontrial{i,f};
%           nose_xcorrs(count,:) = xcf_nose{i,f};
          count = 1+count;
          end
      end
      end
   
   
   
 save(['C:\Users\Nitz_Lab\David\crosscor\' 'xcf_parietal.mat'])
           