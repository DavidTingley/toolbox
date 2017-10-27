function [] = wrapAssemblies_w_theta_phoenix(COND)% try
    % load spikes and stuff
    xml = LoadParameters;
     load([xml.FileName '.behavior.mat'])
    load([xml.FileName '.sessionInfo.mat'])
    lfp = bz_GetLFP(sessionInfo.thetaChans(end));%,'intervals',[behavior.timestamps(1) behavior.timestamps(end)]);
    [b a] = butter(4,[6/(lfp.samplingRate/2) 10/(lfp.samplingRate/2)],'bandpass');
    phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data(:,1)))));
    
    spikes = bz_GetSpikes;
    
    pairs =[];
    for i=1:length((spikes.times))
       for j = i:length((spikes.times)) 
          if i ~= j & spikes.shankID(i) ~= spikes.shankID(j)
              if strcmp(spikes.region{i},'ls') && strcmp(spikes.region{j},'hpc')
                
                   pairs = [pairs; i j]; 
                   
              elseif strcmp(spikes.region{i},'hpc') && strcmp(spikes.region{j},'ls')
                   
                    pairs = [pairs; i j]; 
                    
              end
          end
       end
    end
     if ~exist(['/ifs/data/buzsakilab/results/' spikes.sessionName '_w_positionFits_condition_'  num2str(COND) '.mat'])
    % generate continuous binned spike trains
        for c = 1:length(unique(behavior.events.trialConditions))
            spktrains{c} = [];
            phasetrains{c} =[];
            coords{c} = [];
            velocities{c} = [];
            phasetrains_sin{c} = [];
            phasetrains_cos{c} =[];
        end
        
        for c = 1:size(behavior.events.trialIntervals,1)
              trialLength = ceil(1000*(behavior.events.trialIntervals(c,2)-...
                  behavior.events.trialIntervals(c,1)));
              
              train = single(zeros(length(spikes.times),trialLength));

              for s = 1:length(spikes.times)
                    f = find(spikes.times{s}>behavior.events.trialIntervals(c,1));
                    ff = find(spikes.times{s}<behavior.events.trialIntervals(c,2));
                    fff = intersect(f,ff);
                    times = ceil(1000*(spikes.times{s}(fff)-behavior.events.trialIntervals(c,1)));
                    train(s,times) = 1;
              end
              spktrains{behavior.events.trialConditions(c)} = ...
                  [spktrains{behavior.events.trialConditions(c)},train];
              [a start] =min(abs(lfp.timestamps-behavior.events.trialIntervals(c,1)));
              [a stop] =min(abs(lfp.timestamps-behavior.events.trialIntervals(c,2)));
              p = makeLength(phases(start:stop),length(train));
	      ps = makeLength(sin(phases(start:stop)),length(train));
              pc = makeLength(cos(phases(start:stop)),length(train));


              phasetrains{behavior.events.trialConditions(c)} = ...
                  [phasetrains{behavior.events.trialConditions(c)},p];
              phasetrains_sin{behavior.events.trialConditions(c)} = ...
                  [phasetrains_sin{behavior.events.trialConditions(c)},ps];
              phasetrains_cos{behavior.events.trialConditions(c)} = ...
                  [phasetrains_cos{behavior.events.trialConditions(c)},pc];

              cc = makeLength(behavior.events.trials{c}.mapping,length(train));
              coords{behavior.events.trialConditions(c)} = ...
                  [coords{behavior.events.trialConditions(c)},cc];
              
              vel = abs(diff(behavior.events.trials{c}.x)) + abs(diff(behavior.events.trials{c}.y));
              v = makeLength(vel,length(train));
              velocities{behavior.events.trialConditions(c)} = ...
                  [velocities{behavior.events.trialConditions(c)},v];
            % call GLM here
        end
        for c = COND%length(unique(behavior.events.trialConditions)):-1:1
            extraPredictors = [phasetrains_sin{c};phasetrains_cos{c};phasetrains{c};coords{c};velocities{c}];
            [peerPredictionFits] =cell_assembly_w_theta_crossVal((spktrains{c}),40,extraPredictors,pairs);
            save(['/ifs/data/buzsakilab/results/' spikes.sessionName '_w_positionFits_condition_'  num2str(COND) '.mat'],'peerPredictionFits','-v7.3');
        end
end
%     elseif exist('assembliesCrossRegionData.mat') && ~exist('assemblies.mat')
%         pairs_all = pairs; clear pairs
%         load('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
%         if exist('devControl')
%             for i=1:length(pairs)
%                f = find(sum(ismember(pairs_all,pairs(1,:))')==2);
%                pairs_all(f,:) = [];
%             end
%             pairs = pairs_all;
%             for c = 1:length(unique(behavior.events.trialConditions))
%             spktrains{c} = [];
%             end
%             for c = 1:size(behavior.events.trialIntervals,1)
%             trialLength = ceil(1000*(behavior.events.trialIntervals(c,2)-...
%             behavior.events.trialIntervals(c,1)));
% 
%             train = single(zeros(length(spikes.times),trialLength));
% 
%             for s = 1:length(spikes.times)
%             f = find(spikes.times{s}>behavior.events.trialIntervals(c,1));
%             ff = find(spikes.times{s}<behavior.events.trialIntervals(c,2));
%             fff = intersect(f,ff);
%             times = ceil(1000*(spikes.times{s}(fff)-behavior.events.trialIntervals(c,1)));
%             train(s,times) = 1;
%             end
%             spktrains{behavior.events.trialConditions(c)} = ...
%             [spktrains{behavior.events.trialConditions(c)},train];
%             % call GLM here
%             end
%             for c = 1:length(unique(behavior.events.trialConditions(c)))
%             [dev{c} devControl{c}] =cell_assembly((spktrains{c}),0:200,pairs);
%             save('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
%             end
%         else
%             pairs = pairs_all;
%             for c = 1:size(behavior.events.trialIntervals,1)
%             trialLength = ceil(1000*(behavior.events.trialIntervals(c,2)-...
%             behavior.events.trialIntervals(c,1)));
% 
%             train = single(zeros(length(spikes.times),trialLength));
% 
%             for s = 1:length(spikes.times)
%             f = find(spikes.times{s}>behavior.events.trialIntervals(c,1));
%             ff = find(spikes.times{s}<behavior.events.trialIntervals(c,2));
%             fff = intersect(f,ff);
%             times = ceil(1000*(spikes.times{s}(fff)-behavior.events.trialIntervals(c,1)));
%             train(s,times) = 1;
%             end
%             spktrains{behavior.events.trialConditions(c)} = ...
%             [spktrains{behavior.events.trialConditions(c)},train];
%             % call GLM here
%             end
%             for c = 1:length(unique(behavior.events.trialConditions(c)))
%             [dev{c} devControl{c}] =cell_assembly((spktrains{c}),0:200,pairs);
%             save('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
%             end
%         end
%     end
    
    % save results
%     save('assembliesCrossRegion.mat','dev*','pairs')
% catch
% end
