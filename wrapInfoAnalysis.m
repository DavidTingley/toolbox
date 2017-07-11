


xml = LoadParameters;

load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
load([xml.FileName '.spikes.cellinfo.mat'])
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
nBins = max(behavior.events.trials{1}.mapping);

olypherInfo.region = spikes.region;
olypherInfo.sessionName = spikes.sessionName;
olypherInfo.UID = spikes.UID;
for cell=1:length(spikes.times)
    olypherInfo.results{cell} = table;
end
    

   % set up phase coding data
[rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
[binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior);
    
for smoothing = 1:round(nBins/2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    for cond = 1:length(unique(behavior.events.trialConditions))
%         figure(cond)
        % smooth data..
        for cell = 1:length(spikes.times)
           for trial = 1:size(binnedPhaseMap{cond},2)
              binnedPhaseMap_smooth{cond}(cell,trial,:) = circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,trial,:)),smoothing,'method','mean'); 
              binnedPhaseMap_smooth{cond}(cell,trial,binnedPhaseMap_smooth{cond}(cell,trial,:)==0)=nan;
              rateMap_smooth{cond}(cell,trial,:) = smooth(squeeze(countMap{cond}(cell,trial,:))',smoothing);
              for iter = 1:10
                  binnedPhaseMap_smooth_shuffle{cond}(cell,trial,iter,:) = ...
                      circ_smoothTS(circshift(squeeze(binnedPhaseMap{cond}(cell,trial,:)),round(nBins*rand)),smoothing,'method','mean'); 
                  
                  binnedPhaseMap_smooth_shuffle{cond}(cell,trial,iter,...
                      binnedPhaseMap_smooth_shuffle{cond}(cell,trial,iter,:)==0)=nan;
                  
                  rateMap_smooth_shuffle{cond}(cell,trial,iter,:) = ...
                      smooth(circshift(squeeze(countMap{cond}(cell,trial,:)),round(nBins*rand)),smoothing);
              end
           end
           range = 0:max(max(squeeze(rateMap_smooth{cond}(cell,:,:))))./63:max(max(squeeze(rateMap_smooth{cond}(cell,:,:))));
           if isempty(range)
               range = [0 1];
           end
           % discretize actual data
           f = find(rateMap_smooth{cond}(cell,:)==0);
           rateMap_smooth{cond}(cell,f)=nan;
           f = find(binnedPhaseMap_smooth{cond}(cell,:)==0);
           binnedPhaseMap_smooth{cond}(cell,f)=nan;
           rateMap_disc{cond}(cell,:,:) = discretize(rateMap_smooth{cond}(cell,:,:),range);  % discretize both rate/phase to same # of bins...
           phaseMap_disc{cond}(cell,:,:) = discretize(binnedPhaseMap_smooth{cond}(cell,:,:),-pi:.1:pi);% ,-1:.032:1);
           % discretize shuffled data
           f = find(rateMap_smooth_shuffle{cond}(cell,:)==0);
           rateMap_smooth_shuffle{cond}(cell,f)=nan;
           f = find(binnedPhaseMap_smooth_shuffle{cond}(cell,:)==0);
           binnedPhaseMap_smooth_shuffle{cond}(cell,f)=nan;
           for iter = 1:10
           rateMap_disc_shuffle{cond}(cell,:,iter,:) = discretize(rateMap_smooth_shuffle{cond}(cell,:,iter,:),range);  % discretize both rate/phase to same # of bins...
           phaseMap_disc_shuffle{cond}(cell,:,iter,:) = discretize(binnedPhaseMap_smooth_shuffle{cond}(cell,:,iter,:),-pi:.1:pi);% ,-1:.032:1);
           end
        end
%         phaseMap_disc{cond}(isnan(phaseMap_disc{cond}))=0;    
        
        % run info analysis
        [track_info_rate,pos_info_val_rate] = Info_Analysis(rateMap_disc{cond},1,0);  
        [track_info_phase,pos_info_val_phase] = Info_Analysis(phaseMap_disc{cond},1,0);  
        for iter =1:10
            [track_info_rate_shuffle{iter},pos_info_val_rate_shuffle{iter}] = ...
                Info_Analysis(squeeze(rateMap_disc_shuffle{cond}(:,:,iter,:)),1,0);  
            [track_info_phase_shuffle{iter},pos_info_val_phase_shuffle{iter}] = ...
                Info_Analysis(squeeze(phaseMap_disc_shuffle{cond}(:,:,iter,:)),1,0);  
        end
        
        % compile data
        for cell = 1:length(spikes.times)
            struct.phaseInfoScores = squeeze(sum(pos_info_val_phase(cell,:,:)))';
            struct.phaseTotalInfo = sum(squeeze(sum(pos_info_val_phase(cell,:,:))));
            struct.phasePeakInfo = max(squeeze(sum(pos_info_val_phase(cell,:,:))));
            
            struct.rateInfoScores = squeeze(sum(pos_info_val_rate(cell,:,:)))';
            struct.rateTotalInfo = sum(squeeze(sum(pos_info_val_rate(cell,:,:))));
            struct.ratePeakInfo = max(squeeze(sum(pos_info_val_rate(cell,:,:))));
            
            for iter = 1:10
               struct.rateTotalInfo_shuffle(iter) = sum(squeeze(sum(pos_info_val_rate_shuffle{iter}(cell,:,:))));
               struct.ratePeakInfo_shuffle(iter) = max(squeeze(sum(pos_info_val_rate_shuffle{iter}(cell,:,:))));
               struct.phaseTotalInfo_shuffle(iter) = sum(squeeze(sum(pos_info_val_phase_shuffle{iter}(cell,:,:))));
               struct.phasePeakInfo_shuffle(iter) = max(squeeze(sum(pos_info_val_phase_shuffle{iter}(cell,:,:))));
            end
            struct.smoothing = smoothing;
            struct.condition = cond;
            olypherInfo.results{cell} = [olypherInfo.results{cell}; struct2table(struct)];
%             if cell == 80 & cond == 5
%                 subplot(2,2,1);
%                 imagesc(squeeze(rateMap_disc{cond}(cell,:,:)));
%                 subplot(2,2,2);
%                 imagesc(squeeze(phaseMap_disc{cond}(cell,:,:)));
%                 subplot(2,2,4)
%                 scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end)+2*pi,'.k');
%                 subplot(2,2,3);
%                 rows = find(olypherInfo.results{cell}.condition==cond);
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.ratePeakInfo(rows),'r')
%                 hold on
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.phasePeakInfo(rows),'g')
%                 hold off
%                 title([cell cond])
%                 pause(.1)
%             end
            
        end
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
    end
    olypherInfo.dateRun = date;  % this can take a very long time so lets save each loop...
    save([xml.FileName '.olypherInfo.cellinfo.mat'],'olypherInfo')
end

