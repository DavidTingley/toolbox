


xml = LoadParameters;

load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
load([xml.FileName '.spikes.cellinfo.mat'])
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
bins = max(behavior.events.trials{1}.mapping);

olypherInfo.region = spikes.region;
olypherInfo.sessionName = spikes.sessionName;
olypherInfo.UID = spikes.UID;
for cell=1:length(spikes.times)
    olypherInfo.results{cell} = table;
end
    

   % set up phase coding data
[rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
[binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior);
    
for smoothing = 1:round(bins/2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    
    for cond = 1:length(unique(behavior.events.trialConditions))
        % smooth data..
        for cell = 1:length(spikes.times)
           for trial = 1:size(binnedPhaseMap{cond},2)
              binnedPhaseMap_smooth{cond}(cell,trial,:) = circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,trial,:)),smoothing,'method','mean'); 
              rateMap_smooth{cond}(cell,trial,:) = smooth(countMap{cond}(cell,trial,:),smoothing);
           end
           rateMap_disc{cond}(cell,:,:) = discretize(rateMap_smooth{cond}(cell,:,:),50);  % discretize both rate/phase to same # of bins...
           phaseMap_disc{cond}(cell,:,:) = discretize(binnedPhaseMap_smooth{cond}(cell,:,:),50);
           
        end
        phaseMap_disc{cond}(isnan(phaseMap_disc{cond}))=0;    
        
        % run info analysis
        [track_info_rate,pos_info_val_rate] = Info_Analysis(rateMap_disc{cond},1,0);  
        [track_info_phase,pos_info_val_phase] = Info_Analysis(phaseMap_disc{cond},1,0);  
        
        
        % compile data
        for cell = 1:length(spikes.times)
            struct.phaseInfoScores = squeeze(sum(pos_info_val_phase(cell,:,:)))';
            struct.phaseTotalInfo = sum(squeeze(sum(pos_info_val_rate(cell,:,:))));
            
            struct.rateInfoScores = squeeze(sum(pos_info_val(cell,:,:)))';
            struct.rateTotalInfo = sum(squeeze(sum(pos_info_val(cell,:,:))));
            
            
            struct.smoothing = smoothing;
            struct.condition = cond;
           olypherInfo.results{cell} = [olypherInfo.results{cell}; struct2table(struct)];
        end
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
    end
    olypherInfo.dateRun = date;  % this can take a very long time so lets save each loop...
    save([xml.FileName '.olypherInfo.cellinfo.mat'],'olypherInfo')
end

