function [phaseLockAll] = lfp_spike_phaselock(lfpLocation,spikeLocation,frequencyRange,sampleRate)

% INPUT
% 
%     lfpLocation - folder and filename for the local field potential data
%                   (ex. /home/data/lfp/recording1.mat)
%                   format should be N x T where N is the number of channels
%                   and T is time
%     spikeLocation - folder and filename for the single unit spike times
%                     (ex. /home/data/spikeTrains/recording1.mat)
%                      format should be a cell array {N x 1} where N is the
%                      number of cells and the elements in each cell are spike 
%                      times for that particular cell
%     frequencyRange - (M x 2) range of frequencies(Hz) to analyze 
%                      (ex. [4:100; 5:101]')
%     sampleRate - integer value that represents the rate at which LFP is sampled
%                  (default is 1000 Hz)
% 
% OUTPUT
%     phaseLockAll - 4-d matrix (noChannels x noCells x
%                    size(frequencyRange x 629) where 0:629 corresponds to
%                    0:2pi or 0:360 degrees


if nargin < 4
    sampleRate = 1000;
end

if frequencyRange(1,1) < 4
    error('Frequency range is too low');
elseif frequencyRange(1,1)<10
   filtWinSize = 3;
elseif frequencyRange(1,1)>10
    filtWinSize = 9;
end

if ~ismatrix(spikeLocation)
    spikeTrains = load(spikeLocation);
elseif ismatrix(spikeLocation)
    spikeTrains = spikeLocation;
end
if ~ismatrix(lfpLocation)
    lfp = load(lfpLocation);
elseif ismatrix(lfpLocation)
    lfp = lfpLocation;
end



noCells = length(spikeTrains);
noChannels = size(lfp,1);



for channel = 1:noChannels
    for cell = 1:noCells
        for range = 1:size(frequencyRange,1)
            
            
             [b a] = butter(filtWinSize ,[frequencyRange(range,1)/round(sampleRate/2) ...
                    frequencyRange(range,2)/round(sampleRate/2)],'bandpass');
            
             lfpFiltered = filtfilt(b,a,lfp(channel,:));
             phases = angle(hilbert(lfpFiltered));
             phaseRateWholeRec = zeros(629,1);
              
             for spike = 1:length(spikeTrains{cell})
             
               if length(phases) >= ceil((spikeTrains{cell}(spike)+.0000000001)*1000)
                   phaseRateWholeRec(315+round(phases(ceil((spikeTrains{cell}(spike)+.0000000001)...
                       *1000))*100)) = phaseRateWholeRec(315+round(phases(ceil(...
                       (spikeTrains{cell}(spike)+.0000000001)*1000))*100))+1; 
               elseif length(phases) < ceil(spikeTrains{cell}(spike)*1000)
                   phaseRateWholeRec(315+round(phases(end)*100)) = ...
                        phaseRateWholeRec(315+round(phases(end)*100))+1; 
               end
             
             end
           phaseLockAll(channel,cell,range,:) = phaseRateWholeRec;
        end
    end
end



end