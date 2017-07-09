function [binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior)
% USAGE
% [binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%

nBins = max(behavior.events.trials{1}.mapping);
nCells = length(phaseMap{1});

for cond = 1:length(rateMap)
    
binnedPhaseMap{cond} = nan(size(rateMap{cond},1),size(rateMap{cond},2),size(rateMap{cond},3));

for cell = 1:nCells
    for trial = 1:sum(behavior.events.trialConditions==cond)
        if ~isempty(phaseMap{cond}{cell})
        f = find(phaseMap{cond}{cell}(:,2)==trial);
        if ~isempty(f)
        for bin = 1:nBins
            ff = find(phaseMap{cond}{cell}(f,1)==bin);
             if ~isempty(ff)
                binnedPhaseMap{cond}(cell,trial,bin)=circ_mean(phaseMap{cond}{cell}(f(ff),1));
             end
        end
        end
    end
    end
end
end



end