function [noiseCorrSuccess,noiseCorrFailure,...
          noiseCorrSuccessControl,noiseCorrFailureControl] ...
        = noiseCorr(trials_all_success,trials_all_failure,controlIter,...
          fileName) 
% USAGE
%  [noiseCorrSuccess,noiseCorrFailure,...
%   noiseCorrSuccessControl,noiseCorrFailureControl] ...
%   = noiseCorr(trials_all_success,trials_all_failure,controlIter, fileName)
% 
% INPUT
%         trials_all_success - (N x 1) cell array where each element is a 
%                              matrix (M x D) with M rows corresponding to 
%                              the number of trials and D columns 
%                              corresponding to the number of bins per trial
%
%         trials_all_failure - same as 'trials_all_success' for failure trials
%
%         controlIter - number of iterations to run a boot-strap control
%
%         fileName - optional argument of location and name of 
%                    file you would like results saved to
%                   Example: "C:\Data\noise_correlation_results.mat"
%
%
% OUTPUT
%         noiseCorr* - (N x N) matrices of noise correlation values for
%                      each cell pair
%
%         noiseCorr***Control - (controlIter x N x N) matrices of control noise
%                               correlations where all trials are reorganized 
%                               randomly on each iteration 
% 
% this function calculates the noise correlations for cell pairs across a set of trials,
% relative to time or space. 
%        
% David Tingley, 2015

totalCellNumber = length(trials_all_success);
disp(totalCellNumber)
if totalCellNumber > 1
    if nargin < 3
        controlIter = 50;
    end

    for numCell = 1:totalCellNumber
        success_mean = mean(trials_all_success{numCell});
        failure_mean = mean(trials_all_failure{numCell});

        numTrialsSuccess = size(trials_all_success{numCell},1);
        numTrialsFailure = size(trials_all_failure{numCell},1);
        numBins = size(trials_all_success{numCell},2);

        % subtract mean rate from individual trials and reshape into (Z x 1)
        % vector 
        successDifferences{numCell} = reshape((trials_all_success{numCell}(:,:)-...
            repmat(success_mean,size(trials_all_success{numCell},1),1))',...
            numTrialsSuccess*numBins,1);
        failureDifferences{numCell} = reshape((trials_all_failure{numCell}(:,:)-...
            repmat(failure_mean,size(trials_all_failure{numCell},1),1))',...
            numTrialsFailure*numBins,1);       


        %% Control by randomizing trial order and re-running
            for p = 1:controlIter
                % boot-strap control for success trials
                r = randperm(size(trials_all_success{numCell},1));

            successDifferencesControl{p,numCell} = reshape((trials_all_success{numCell}(r,:)-...
            repmat(success_mean,size(trials_all_success{numCell},1),1))',...
            numTrialsSuccess*numBins,1);
                % boot-strap control for failure trials
                r = randperm(size(trials_all_failure{numCell},1));

            failureDifferencesControl{p,numCell} = reshape((trials_all_failure{numCell}(r,:)-...
            repmat(failure_mean,size(trials_all_failure{numCell},1),1))',...
            numTrialsFailure*numBins,1);       
            end
    end

        for i =  1:totalCellNumber
        for j =  1:totalCellNumber
            if j ~= i 
                if size(successDifferences{i})==size(successDifferences{j})
                    if size(failureDifferences{i})==size(failureDifferences{j})
                % correlate variation from mean("noise") for each cell pair
                % across all trials
                tempSuccess  = corrcoef(successDifferences{i},successDifferences{j});
                tempFailure = corrcoef(failureDifferences{i},failureDifferences{j});
                noiseCorrSuccess(i,j)  = tempSuccess(2);
                noiseCorrFailure(i,j)  = tempFailure(2);

                for p = 1:controlIter
                tempSuccessControl = corrcoef(successDifferencesControl{p,i},...
                    successDifferencesControl{p,j});
                tempFailureControl = corrcoef(failureDifferencesControl{p,i},...
                    failureDifferencesControl{p,j});

                noiseCorrSuccessControl(p,i,j)  = tempSuccessControl(2);
                noiseCorrFailureControl(p,i,j)  = tempFailureControl(2);
                end

                    end
                end
            end
        end
        end

        if nargin > 3
            save([fileName],'-v7.3')
        end
end

 return
 
 
 
 
 
