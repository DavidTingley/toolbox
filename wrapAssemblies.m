try
    % load spikes and stuff
    xml = LoadParameters;
     load([xml.FileName '.behavior.mat'])
    load([xml.FileName '.sessionInfo.mat'])
    
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
    if ~exist('assembliesCrossRegionData.mat') && ~exist('assemblies.mat')
    % generate continuous binned spike trains
        for c = 1:length(unique(behavior.events.trialConditions))
            spktrains{c} = [];
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
            % call GLM here
        end
        for c = 1:length(unique(behavior.events.trialConditions))
            [dev{c} devControl{c}] =cell_assembly((spktrains{c}),0:200,pairs);
            save('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
        end
    elseif exist('assembliesCrossRegionData.mat') && ~exist('assemblies.mat')
        pairs_all = pairs; clear pairs
        load('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
        if exist('devControl')
            for i=1:length(pairs)
               f = find(sum(ismember(pairs_all,pairs(1,:))')==2);
               pairs_all(f,:) = [];
            end
            pairs = pairs_all;
            for c = 1:length(unique(behavior.events.trialConditions))
            spktrains{c} = [];
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
            % call GLM here
            end
            for c = 1:length(unique(behavior.events.trialConditions(c)))
            [dev{c} devControl{c}] =cell_assembly((spktrains{c}),0:200,pairs);
            save('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
            end
        else
            pairs = pairs_all;
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
            % call GLM here
            end
            for c = 1:length(unique(behavior.events.trialConditions(c)))
            [dev{c} devControl{c}] =cell_assembly((spktrains{c}),0:200,pairs);
            save('assembliesCrossRegionData.mat','spktrains','pairs','dev*');
            end
        end
    end
    
    % save results
%     save('assembliesCrossRegion.mat','dev*','pairs')
catch
end