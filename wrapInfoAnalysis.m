


% try
xml = LoadParameters;
% if exist([xml.FileName '.olypherInfo_w_disc.cellinfo.mat'])
load([xml.FileName '.firingMaps.cellinfo.mat'])
load([xml.FileName '.phaseMaps.cellinfo.mat'])
load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
load([xml.FileName '.spikes.cellinfo.mat'])
% lfp = bz_GetLFP(sessionInfo.thetaChans(2));
nBins = max(behavior.events.trials{1}.mapping);

olypherInfo.region = spikes.region;
olypherInfo.sessionName = spikes.sessionName;
olypherInfo.UID = spikes.UID;
for cell=1:length(spikes.times)
    olypherInfo.results{cell} = table;
end
    

%    % set up phase coding data
% [firingMaps.rateMaps firingMaps.countMaps occuMap firingMaps.phaseMaps] = bz_firingMap1D(spikes.times,behavior,lfp,4);
[binnedfiringMaps.phaseMaps] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
    for discBins = [1:4]%1:8 10 15 40 60 100]
for smoothing = [1:10 20 50]%1:round(nBins/2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    for cond = 1:length(unique(behavior.events.trialConditions))
%         figure(cond)
        % smooth data..
        if sum(behavior.events.trialConditions==cond)>2
        for cell = 1:length(spikes.times)
           for trial = 1:size(binnedfiringMaps.phaseMaps{cond},2)
              binnedfiringMaps.phaseMaps_smooth{cond}(cell,trial,:) = circ_smoothTS(squeeze(binnedfiringMaps.phaseMaps{cond}(cell,trial,:)),smoothing,'method','mean','exclude',0); 
              binnedfiringMaps.phaseMaps_smooth{cond}(cell,trial,binnedfiringMaps.phaseMaps_smooth{cond}(cell,trial,:)==0)=nan;
              firingMaps.rateMaps_smooth{cond}(cell,trial,:) = smooth(squeeze(firingMaps.countMaps{cond}(cell,trial,:))',smoothing);
%               for iter = 1:10
%                   binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,trial,iter,:) = ...
%                       circ_smoothTS(circshift(squeeze(binnedfiringMaps.phaseMaps{cond}(cell,trial,:)),round(nBins*rand)),smoothing,'method','mean'); 
%                   
%                   binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,trial,iter,...
%                       binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,trial,iter,:)==0)=nan;
%                   
%                   firingMaps.rateMaps_smooth_shuffle{cond}(cell,trial,iter,:) = ...
%                       smooth(circshift(squeeze(firingMaps.countMaps{cond}(cell,trial,:)),round(nBins*rand)),smoothing);
%               end
           end
           range = 0:max(max(squeeze(firingMaps.rateMaps_smooth{cond}(cell,:,:))))./discBins:max(max(squeeze(firingMaps.rateMaps_smooth{cond}(cell,:,:))));
           if isempty(range)
               range = [0 1];
           end
           % discretize actual data
           f = find(firingMaps.rateMaps_smooth{cond}(cell,:)==0);
           firingMaps.rateMaps_smooth{cond}(cell,f)=nan;
           f = find(binnedfiringMaps.phaseMaps_smooth{cond}(cell,:)==0);
           binnedfiringMaps.phaseMaps_smooth{cond}(cell,f)=nan;
           firingMaps.rateMaps_disc{cond}(cell,:,:) = discretize(firingMaps.rateMaps_smooth{cond}(cell,:,:),range);  % discretize both rate/phase to same # of bins...
           firingMaps.phaseMaps_disc{cond}(cell,:,:) = discretize(binnedfiringMaps.phaseMaps_smooth{cond}(cell,:,:),-pi:(2*pi)./discBins:pi);% ,-1:.032:1);
           
           % discretize shuffled data
%            f = find(firingMaps.rateMaps_smooth_shuffle{cond}(cell,:)==0);
%            firingMaps.rateMaps_smooth_shuffle{cond}(cell,f)=nan;
%            f = find(binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,:)==0);
%            binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,f)=nan;
%            for iter = 1:10
%            firingMaps.rateMaps_disc_shuffle{cond}(cell,:,iter,:) = discretize(firingMaps.rateMaps_smooth_shuffle{cond}(cell,:,iter,:),range);  % discretize both rate/phase to same # of bins...
%            firingMaps.phaseMaps_disc_shuffle{cond}(cell,:,iter,:) = discretize(binnedfiringMaps.phaseMaps_smooth_shuffle{cond}(cell,:,iter,:),-pi:.1:pi);% ,-1:.032:1);
%            end
        end
%         firingMaps.phaseMaps_disc{cond}(isnan(firingMaps.phaseMaps_disc{cond}))=0;    
        
        % run info analysis
        [track_info_rate,pos_info_val_rate] = Info_Analysis(firingMaps.rateMaps_disc{cond},1,0);  
        [track_info_phase,pos_info_val_phase] = Info_Analysis(firingMaps.phaseMaps_disc{cond},1,0);  
%         for iter =1:10
%             [track_info_rate_shuffle{iter},pos_info_val_rate_shuffle{iter}] = ...
%                 Info_Analysis(squeeze(firingMaps.rateMaps_disc_shuffle{cond}(:,:,iter,:)),1,0);  
%             [track_info_phase_shuffle{iter},pos_info_val_phase_shuffle{iter}] = ...
%                 Info_Analysis(squeeze(firingMaps.phaseMaps_disc_shuffle{cond}(:,:,iter,:)),1,0);  
%         end
        
        % compile data
        for cell = 1:length(spikes.times)
            struct.phaseInfoScores = squeeze(sum(pos_info_val_phase(cell,:,:)))';
            struct.phaseTotalInfo = sum(squeeze(sum(pos_info_val_phase(cell,:,:))));
            struct.phasePeakInfo = max(squeeze(sum(pos_info_val_phase(cell,:,:))));
%             phase = max(squeeze(sum(pos_info_val_phase(cell,:,:))))
            struct.rateInfoScores = squeeze(sum(pos_info_val_rate(cell,:,:)))';
            struct.rateTotalInfo = sum(squeeze(sum(pos_info_val_rate(cell,:,:))));
            struct.ratePeakInfo = max(squeeze(sum(pos_info_val_rate(cell,:,:))));
%             rate = max(squeeze(sum(pos_info_val_rate(cell,:,:))))
%             for iter = 1:10
%                struct.rateTotalInfo_shuffle(iter) = sum(squeeze(sum(pos_info_val_rate_shuffle{iter}(cell,:,:))));
%                struct.ratePeakInfo_shuffle(iter) = max(squeeze(sum(pos_info_val_rate_shuffle{iter}(cell,:,:))));
%                struct.phaseTotalInfo_shuffle(iter) = sum(squeeze(sum(pos_info_val_phase_shuffle{iter}(cell,:,:))));
%                struct.phasePeakInfo_shuffle(iter) = max(squeeze(sum(pos_info_val_phase_shuffle{iter}(cell,:,:))));
%             end
            struct.smoothing = smoothing;
            struct.discBins = discBins;
            struct.condition = cond;
            olypherInfo.results{cell} = [olypherInfo.results{cell}; struct2table(struct)];
%             if cell == 80 && cond == 5
%                 rows = find(olypherInfo.results{cell}.condition==cond);
%                 cols = find(olypherInfo.results{cell}.discBins == discBins);
%                 rows = intersect(rows,cols);
%                 figure(discBins)
%                 subplot(4,2,1);
% %                 imagesc(squeeze(firingMaps.rateMaps_disc{cond}(cell,:,:)));
%                 imagesc(squeeze(firingMaps.phaseMaps_disc{cond}(cell,:,:)));
%                 subplot(4,2,2);
%                 rows = find(olypherInfo.results{cell}.condition==cond);
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.rateTotalInfo(rows),'r')
%                 hold on
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.phaseTotalInfo(rows),'g')
%                 hold off
%                 subplot(4,2,4)
%                 scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),phaseMaps.phaseMaps{cond}{cell}(:,end)+2*pi,'.k');
%                 subplot(4,2,3);
%                 rows = find(olypherInfo.results{cell}.condition==cond);
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.ratePeakInfo(rows),'r')
%                 hold on
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.phasePeakInfo(rows),'g')
%                 hold off
%                 subplot(4,2,5)
%                 imagesc(olypherInfo.results{cell}.phaseInfoScores)
%                 subplot(4,2,6)
%                 imagesc(olypherInfo.results{cell}.rateInfoScores)
%                 
%                 title([cell cond])
%                 pause(.1)
%             end
            
        end
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
        end
    end
    olypherInfo.dateRun = date;  % this can take a very long time so lets save each loop...
    save([xml.FileName '.olypherInfo_w_disc.cellinfo.mat'],'olypherInfo','-v7.3')
% end
    end
% end
% catch
end