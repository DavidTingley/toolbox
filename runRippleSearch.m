% 
% PulsePal
ProgramPulsePalParam(1, 14, 1); % Sets output channel 1 to use custom train 1
ProgramPulsePalParam(2, 14, 1); % Sets output channel 1 to use custom train 1
ProgramPulsePalParam(3, 14, 2); % Sets output channel 1 to use custom train 1
% %             ProgramPulsePalParam(4, 14, 1); % Sets output channel 1 to use custom train 1
% 
% 
% % SineWaveVoltages = sin((2*pi/100)*(1:1000));
rateRange = linspace(0.02,2,6);
% rateRange = 2.5;rateRange(randperm(length(rateRange)));
% 
% amplitudes = linspace(1,5.9,6); %1.1 1.4 1.64
% amplitudes = amplitudes(randperm(length(amplitudes)));
% 
% amplitudes = 4;5.9;
% 
% window = 200; %1:25:300 divide by 10 to get stim duration
% amplitude = 4;
amplitudes = linspace(1.5,4,8);

% for iter = 1:100
    tic
for a = 1:length(amplitudes)
    amplitude = amplitudes(a);
    waveform = [zeros(50,1); amplitude*ones(window,1); zeros(50,1)]';
    
%     waveform = amplitude * minmax(1:window);
    waveform2 = waveform ./ 10;
    NOTE = 'analogin pulses are 1/10th actual voltage';
    %             waveform =  1.+.5*minmax_norm(1:800);


    % NoiseWaveVoltages = randn(1,1000)*2;

    % For the next six lines, note the parameter codes for ProgramPulsePalParam.

    SendCustomWaveform(1, 0.0001, waveform); % Uploads sine waveform. Samples are played at 10khz.
    SendCustomWaveform(2, 0.0001, waveform2); % Uploads sine waveform. Samples are played at 10khz.


for r = length(rateRange):-1:1
    disp(['running with rate: ' num2str(rateRange(r))])
    nPulses = ceil(rateRange(r)*60);
    pattern{r} = randi([20,round(2/rateRange(r)*1000)],ceil(nPulses),1); % create IRI in milliseconds
        
    while nPulses > 0
    %             waveform = amplitude * minmax(1:window*2);
%             waveform = amplitude * minmax_norm(Smooth([zeros(499,1); 1; zeros(500,1)],window))';
          
            TriggerPulsePal('0111'); % Soft-triggers channels 1, 3 and 4.
            disp(['delay: ' num2str(pattern{r}(ceil(nPulses))/1000) ', ' num2str(ceil(nPulses)) ' pulses left...'])
            pause(pattern{r}(ceil(nPulses))/1000)
            nPulses = nPulses - 1;
    end
    pause(60*randi(6))
end
% end

toc
disp('done')

save(['pupilStim_' datestr(now,30) '.mat'])
end
% 2. (Example requiring the longer 5,000 pulse sequences on Pulse Pal 2)

% Attach a speaker to output channel 1, to hear the waveform.

% Wave = load ('PulsePalWave.mat'); Wave = Wave.PulsePalWave; % Loads a complex voltage waveform
% 
% ProgramPulsePalParam(1, 'CustomTrainID', 1); % Sets output channel 1 to use custom train 1
% 
% SendCustomWaveform(1, 0.0001, Wave); % Sends the waveform to Pulse Pal, sampled at 10kHz
% 
% TriggerPulsePal(1); % Triggers channel 1. Make sure your speaker volume is audible.

% EndPulsePal