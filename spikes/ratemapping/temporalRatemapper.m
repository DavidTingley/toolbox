function temporalRatemapper(sigfileList,eventTS,trialTypes,binSize,dataLocation,ratemapFolder,ratemapName)
% INPUTS
%     sigfileList - cell array of ascii filenames that contain
%                   spike time information (spike timestamps)
%     eventTS - matrix (M x N) where M rows signify each trial that was 
%               recorded and N columns represent timestamps signifying 
%               particular events within each trial (i.e. stimuli 
%               presentation or reward)
%     trialTypes - vector (M x 1) where each element depicts which type of 
%                  trial (i.e. success, failure, or quinine)
%     binSize - integer value that signifies the desired bin size to be 
%               ratemaped (this will be an approximation depending on the
%               mean/median length of time between behavioral events in 
%               'eventTS') 
%     ratemapFolder - folder location where you would like final ratemapping 
%                     output to be saved
%     dataLocation - folder location where ascii sigfiles are saved
%     ratemapName - name of ratemapping file you would like to be saved
% SAVED OUTPUTS
%



numEvents = size(eventTS,2);


for p = 1:numEvents-1
   meanPhaseLength(p) = mean(eventTS(:,p+1)-eventTS(:,p))*1000;
   medianPhaseLength(p) = median(eventTS(:,p+1)-eventTS(:,p))*1000;
   
   numBins(p) = round(meanPhaseLength(p)/binSize);
end


for type = 1:max(trialTypes)
    trials = find(trialTypes==type);
    
    [rate_smooth{type},rate_smooth_mean{type},rate_mean{type}] =...
    ratemap(sigfileList,dataLocation,eventTS(trials,:),numBins)
end

save([ratemapFolder '/' ratemapName '.mat'])

return

function [rate_smooth,rate_smooth_mean,rate_mean] =...
    ratemap(sigfileList,dataLocation,eventTS,numBins)

for numCell = 1:length(sigfileList)
%     cellTS = load([dataLocation '/' sigfileList{numCell}]);
    cellTS = sigfileList{numCell}(:,1);
    for t = 1:size(eventTS,1)
      c = 1;
      
    for p = 1:length(numBins)
           trial_timebin = (eventTS(t,p+1)-eventTS(t,p))/numBins(p);
           cc = 1;
       for b = 1:numBins(p)   
       numSpikes =  find(cellTS > eventTS(t,p)+cc*trial_timebin);
       numSpikes2 = find(cellTS <= eventTS(t,p)+cc*trial_timebin + (trial_timebin - .000001));
       numSpikes = intersect(numSpikes,numSpikes2);
       rate(numCell,t,c) = length(numSpikes)*(1/trial_timebin);
       c = 1+c;
       cc=cc+1;
       end

    end
    
    z_peth=0:1:4;
    y3_peth=gaussmf(z_peth,[1 3]);
    y3_peth=y3_peth./sum(y3_peth);

    rate_smooth(numCell,t,:) = filtfilt(y3_peth,1,rate(numCell,t,:));

    end
    rate_smooth_mean(numCell,:) = squeeze(mean(rate_smooth(numCell,:,:),2));
    rate_mean(numCell,:) = squeeze(mean(rate(numCell,:,:),2));
end

return
