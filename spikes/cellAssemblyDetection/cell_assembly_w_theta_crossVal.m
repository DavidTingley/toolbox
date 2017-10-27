function [peerPredictionFits] =... %stats devControl statsControl] = ...
    cell_assembly_w_theta_crossVal(spikeTimes,winRange,extraPredictors,pairsToRun)

% INPUT
% 
% 
%     spikeTimes - (M x N x D) matrix of spike times. M is the cell number, 
%                  N is the trial number, and D is length of a trial in
%                  milliseconds
%     winRange - vector that represents the range of temporal smoothing 
%                 windows over which to run the assembly analysis 
%                 (default is 0:500 ms)
%     pairsToRun - (P x 2) matrix of cell pairs that will be analyzed 
%                  (default analyzes all cell pairs within spikeTimes matrix)

% OUTPUT
%      dev - (W x M x N) matrix of deviance from solution vectors for each 
%              cell for each trial at each temporally smoothed window size
%      stats - cell array of statistics returns by glmfit.m at each 
%                temporally smoothed window size
%      *Control - same variables as above with trial order randomly shuffled





if nargin < 1
    error('Must have at least one argument: cell_assembly(spikeTimes)')
elseif nargin == 1
    numCells = size(spikeTimes,1);
    winRange = 0:500;
    pairsToRun = [];
    for n = 1:numCells
        for nn = n:numCells
            if n ~= nn
        pairsToRun = [pairsToRun; n,nn];
            end
        end
    end
elseif nargin == 2
    numCells = size(spikeTimes,1);
    pairsToRun = [];
    for n = 1:numCells
        for nn = n:numCells
            if n ~= nn
        pairsToRun = [pairsToRun; n,nn];
            end
        end
    end
end


numTrials = size(spikeTimes,2);
if length(size(spikeTimes))==2
   numTrials = 1;
   spikes(:,1,:) = spikeTimes; 
   spikeTimes = spikes;
end
% EDIT THIS to take into acount matrices with a single trial (M x D) rather
% than (M x N x D)

% NOTE: Indeces in all outputs match the intervals of winRange, which is not
%       necessarily intervals of 1.
c = 1;

for win = winRange
%     tic
    %              stats{c}   devControl(c,:,:) statsControl{c}
       [peerPredictionFits(c)] = ...
           parGLMRun(win,spikeTimes,extraPredictors,pairsToRun,numTrials);
       c = c+1
       win
%        toc
end
end

%               stats devControl statsControl
function [peerPredictionFits] = ...
    parGLMRun(win,spikeTimes,extraPredictors,pairsToRun,numTrials)
warning off;
peerPredictionFits.results = table;
pred_last = 0;
for pair = 1:size(pairsToRun(:,1),1)
%     tic
pred = pairsToRun(pair,1);
act = pairsToRun(pair,2);
if sum(spikeTimes(act,:,:))+sum(spikeTimes(pred,:,:)) > 100
   %% create temporally smoothed predictor spike trains here
if pred_last ~= pred
    for temp = 1:size(spikeTimes,3)
        if temp>win & temp+win<size(spikeTimes,3)
            if win == 0 
            smoothedTrains(:,temp)= spikeTimes(pred,:,temp);
            elseif win > 0
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-win:temp+win))/(win*2),3);
            end
        end
        if temp<=win
            j = temp-1;
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-j:temp+win))/(length(1:temp+win)),3);
        end
        if temp+win>=size(spikeTimes,3)
            j = size(spikeTimes,3)-temp;
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-win:temp+j))/(length(temp-win:size(spikeTimes,3))),3);
        end
%         smoothedTrains(:,temp) = smooth(spikeTimes(pred,:,temp),win);
    end
end
         %% Run GLM on individual trials/instances
        actual = [];
        predictor =[];
        for trial = 1:numTrials
            actual = [actual;squeeze(double(squeeze(spikeTimes(act,trial,:))))]; % GLM's don't like singles...?
            predictor = [predictor;squeeze(double(smoothedTrains(trial,:)))];
        end
        struct.tau = win;
        
        for iteration = 1:10
            % 60/40 split here
            r = randperm(length(actual));

            actual_test = actual(r(prctile(1:length(actual),60):end));
            predictor_test = predictor(r(prctile(1:length(actual),60):end));
            actual_train = actual(r(1:prctile(1:length(actual),60)));
            predictor_train = predictor(r(1:prctile(1:length(actual),60)));
            extraPredictors_test = extraPredictors(:,r(prctile(1:length(actual),60):end));
            extraPredictors_train = extraPredictors(:,r(1:prctile(1:length(actual),60)));


            [results dev ] = ...
            glmfit([predictor_train;extraPredictors_train]',actual_train,'normal');
            yhat = glmval(results,[predictor_test;extraPredictors_test]','identity');

            rr = randperm(length(actual_train));
            rrr = randperm(length(actual_test));
            [results_shuffle dev_shuffle ] = ...
            glmfit([predictor_train;extraPredictors_train]',actual_train(rr),'normal');
            yhat_shuffle = glmval(results,[predictor_test(rrr);extraPredictors_test(:,rrr)]','identity');
            
            struct.dev = dev;
            struct.dev_shuffle = dev_shuffle;
            struct.iteration=iteration;
            fits = cell2struct({yhat,actual_test,yhat_shuffle,extraPredictors_test(4,:)'},{'yfit','rate','yfit_shuffle','position'},2);
            struct.fits = fits;
            struct.pair = pair;
            struct.ls = pairsToRun(pair,1);
            struct.hpc = pairsToRun(pair,2);
            peerPredictionFits.results = [peerPredictionFits.results;struct2table(struct)];
        end
            
         
     pred_last = pred;
%      toc
 end
end
end
