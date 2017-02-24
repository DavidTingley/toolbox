function [] = lfp_ratemapper_parietal()
% clear all;

rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];

% rstart = 1;

% load('/home/uqdtingl/Dropbox/Nitz lab Shared/TBG_paper/population_rate_lfp_data_f.mat','files')
load('/home/david/Dropbox/Documents/pubs/bf_theta_beta_gammas/file.mat')
load('/home/david/Dropbox/Documents/pubs/bf_theta_beta_gammas/S_markers_times.mat')
load('/zpool/nitzlabdata/F_markers_times.mat')
cc=1;
% 
% for r = rstart:8
% 

    list = dir(['/home/david/Dropbox/LFP_parietal/'  '*.mat']);
 
    
for i = 1:length(list)
   
         

   load(['/home/david/Dropbox/LFP_parietal/' list(i).name],'ad');
%    ad = add;
c=1;
while isempty(ad{c})
    c=1+c;
    if c > 20
       error('dah fu') 
    end
end
ad = ad{c};

   sRec=[];
   for t = 1:length(file)
      sName = [file(t).name(1:7) '_' file(t).name(8:9)];
      lName = list(i).name(1:10);
      if strcmp(sName,lName)
         sRec = t; 
      end
   end


[b a] = butter(3,[4/200 9/200],'bandpass');
temp_wave=[];
temp_wave_control=[];       
pa_splined = [];
% if ~exist(['/zpool/nitzlabdata/LFP_parietal_data_phase_amp/' list(i).name])  
   if ~isempty(sRec)
%    for j = 1:size(ad,1)
       if ~isnan(ad) & ~isempty(ad)


        %% phase/amplitude stuff here
       freq = scal2frq(1:150,'cgau4',1/1000);

        c = cwt(ad,freq,'cgau4');
        S = abs(c.*c); clear c;
        wavelet_all = [(100*S./sum(S(:)))];
        phases = angle(hilbert(filtfilt(b,a,ad)));
        clear S
% successful trials
       
       
                for trial = 1:size(S_markers_times{sRec},1)
                    
                    start = ceil((S_markers_times{sRec}(trial,1))*1000);
                    light = ceil((S_markers_times{sRec}(trial,2))*1000);
                    nosepoke = ceil((S_markers_times{sRec}(trial,3))*1000);
                    platecross = ceil((S_markers_times{sRec}(trial,4))*1000);
                    stop = ceil((S_markers_times{sRec}(trial,5))*1000);
                    trialend = ceil((S_markers_times{sRec}(trial,6))*1000);
                    if trialend < size(wavelet_all,2)
                        
                     %% phase/amplitude stuff here
                    pa{trial} = zeros(trialend-start,3,35);
                    for p = start:trialend
                        for h = 1:35
                           f = find(round((phases(p-250:p+250)+pi)*100/18)==h);
                           
                           pa{trial}(p-start+1,1,h) = mean(mean(wavelet_all(20:35,f)'));
                           pa{trial}(p-start+1,2,h) = mean(mean(wavelet_all(45:65,f)'));
                           pa{trial}(p-start+1,3,h) = mean(mean(wavelet_all(80:150,f)'));
                        end
                    end
                    for tt = 1:35
                    add = 0;
                    len =  [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./119.5:size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./79.5:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./159.5:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./59.5:size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,1,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./119.5:size(pa{trial}(stop+1-start:trialend-start),2))];
                    while length(len) > 540
                        add = add-.01;
                         len =  [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,1,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{trial}(stop+1-start:trialend-start),2))];
                    end
                    while length(len) < 540
                        add = add+.01;

                        len =  [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,1,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{trial}(stop+1-start:trialend-start),2))];
                    end
                    
                    
                    
                        pa_splined(trial,1,:,tt) = ...
                    [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,1,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{trial}(stop+1-start:trialend-start),2))];
                
                    pa_splined(trial,2,:,tt) = ...
                    [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,2,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,2,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,2,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,2,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,2,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{trial}(stop+1-start:trialend-start),2))];
                
                
                    pa_splined(trial,3,:,tt) = ...
                    [spline(1:size(pa{trial}(start+1-start:light+1+1-start),2),pa{trial}(start+1-start:light+1+1-start,3,tt),1:size(pa{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2),pa{trial}(light+1-start:nosepoke+1+1-start,3,tt),1:size(pa{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2),pa{trial}(nosepoke+1-start:platecross+1+1-start,3,tt),1:size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa{trial}(platecross+1-start:stop+1+1-start),2),pa{trial}(platecross+1-start:stop+1+1-start,3,tt),1:size(pa{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa{trial}(stop+1-start:trialend-start),2),pa{trial}(stop+1-start:trialend-start,3,tt),1:size(pa{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{trial}(stop+1-start:trialend-start),2))];
                
                    end
                    %%
  
                    for tt = 1:150
                    temp_wave(trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                     r = randi(10000)-5000;
                    start = ceil((S_markers_times{sRec}(trial,1))*1000)+r;
                    light = ceil((S_markers_times{sRec}(trial,2))*1000)+r;
                    nosepoke = ceil((S_markers_times{sRec}(trial,3))*1000)+r;
                    platecross = ceil((S_markers_times{sRec}(trial,4))*1000)+r;
                    stop = ceil((S_markers_times{sRec}(trial,5))*1000)+r;
                    trialend = ceil((S_markers_times{sRec}(trial,6))*1000)+r;
                    
                    
                      for tt = 1:150
                    temp_wave_control(trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                      end
                    
%                     for tt = 1:150
%                     temp_wave2(trial,tt,:) = ...
%                     [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:120),...
%                     spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:80),...
%                     spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:160),...
%                     spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:60),...
%                     spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:120)];
%                     end
                    end

                end
        
  % failure trials
                  temp_wave_f=[];
                pa_splined_f=[];
                for trial = 1:size(F_markers_times{sRec},1)
                    
                    start = ceil((F_markers_times{sRec}(trial,1))*1000);
                    light = ceil((F_markers_times{sRec}(trial,2))*1000);
                    nosepoke = ceil((F_markers_times{sRec}(trial,3))*1000);
                    platecross = ceil((F_markers_times{sRec}(trial,4))*1000);
                    stop = ceil((F_markers_times{sRec}(trial,5))*1000);
                    trialend = ceil((F_markers_times{sRec}(trial,6))*1000);
                    if trialend < size(wavelet_all,2)
                        if trialend-start > 2200
                        
                     %% phase/amplitude stuff here
                    pa_f{trial} = zeros(trialend-start,3,35);
                    for p = start:trialend
                        for h = 1:35
                           f = find(round((phases(p-250:p+250)+pi)*100/18)==h);
                           
                           pa_f{trial}(p-start+1,1,h) = mean(mean(wavelet_all(20:35,f)'));
                           pa_f{trial}(p-start+1,2,h) = mean(mean(wavelet_all(45:65,f)'));
                           pa_f{trial}(p-start+1,3,h) = mean(mean(wavelet_all(80:150,f)'));
                           if isempty(f)
                               pa_f{trial}(p-start+1,1,h) = nan;
                               pa_f{trial}(p-start+1,2,h) = nan;
                               pa_f{trial}(p-start+1,3,h) = nan;
                           end
                        end
                    end
                    for tt = 1:35
                    add = 0;
                    len =  [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./119.5:size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./79.5:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./159.5:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./59.5:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,1,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./119.5:size(pa_f{trial}(stop+1-start:trialend-start),2))];
                    while length(len) > 540
                        add = add-.01;
                         len =  [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,1,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{trial}(stop+1-start:trialend-start),2))];
                    end
                    while length(len) < 540
                        add = add+.01;

                        len =  [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,1,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{trial}(stop+1-start:trialend-start),2))];
                    end
                    
                    
                    
                        pa_splined_f(trial,1,:,tt) = ...
                    [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,1,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,1,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{trial}(stop+1-start:trialend-start),2))];
                
                    pa_splined_f(trial,2,:,tt) = ...
                    [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,2,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,2,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,2,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,2,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,2,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{trial}(stop+1-start:trialend-start),2))];
                
                
                    pa_splined_f(trial,3,:,tt) = ...
                    [spline(1:size(pa_f{trial}(start+1-start:light+1+1-start),2),pa_f{trial}(start+1-start:light+1+1-start,3,tt),1:size(pa_f{trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{trial}(start+1-start:light+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2),pa_f{trial}(light+1-start:nosepoke+1+1-start,3,tt),1:size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{trial}(light+1-start:nosepoke+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{trial}(nosepoke+1-start:platecross+1+1-start,3,tt),1:size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{trial}(nosepoke+1-start:platecross+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2),pa_f{trial}(platecross+1-start:stop+1+1-start,3,tt),1:size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{trial}(platecross+1-start:stop+1+1-start),2)),...
                    spline(1:size(pa_f{trial}(stop+1-start:trialend-start),2),pa_f{trial}(stop+1-start:trialend-start,3,tt),1:size(pa_f{trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{trial}(stop+1-start:trialend-start),2))];
                
                    end
                    %%
  
                    for tt = 1:150
                    temp_wave_f(trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                    
                      r = randi(5000)-2500;
                    start = ceil((F_markers_times{sRec}(trial,1))*1000)+r;
                    light = ceil((F_markers_times{sRec}(trial,2))*1000)+r;
                    nosepoke = ceil((F_markers_times{sRec}(trial,3))*1000)+r;
                    platecross = ceil((F_markers_times{sRec}(trial,4))*1000)+r;
                    stop = ceil((F_markers_times{sRec}(trial,5))*1000)+r;
                    trialend = ceil((F_markers_times{sRec}(trial,6))*1000)+r;
                    
                        for tt = 1:150
                    temp_wave_control_f(trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                        end
                    
                    
%                     for tt = 1:150
%                     temp_wave2(trial,tt,:) = ...
%                     [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:119.5),...
%                     spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:79.5),...
%                     spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:159.5),...
%                     spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:59.5),...
%                     spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:119.5)];
%                     end
                        end
                    end

                end
           
   

       end
       
%    end
   save(['/zpool/nitzlabdata/LFP_parietal_data_phase_amp/' list(i).name],'wavelet_all','pa*','temp*','-v7.3')
%    end
   
%     end
    
matchedRecs(cc,:) = [i,sRec];
lfp_rec_wave{cc} = (temp_wave);
lfp_rec_wave_control{cc} = (temp_wave_control);
lfp_rec_wave_control_f{cc} = (temp_wave_control_f);
lfp_rec_pa{cc} = (pa_splined);

lfp_rec_wave_f{cc} = (temp_wave_f);

 lfp_rec_pa_f{cc} = (pa_splined_f);
clear pa*
cc=1+cc

end



end

   save(['/zpool/nitzlabdata/LFP_parietal_wavelet_phase_amp_data.mat'],'-v7.3')

% save('/home/uqdtingl/Dropbox/Nitz lab Shared/TBG_paper/lfp_across_trials.mat','-v7.3')

end