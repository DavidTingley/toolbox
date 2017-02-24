function [] = lfp_ratemapper()
% clear all;

rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
dbstop if error

% rstart = 1;

% load('/home/uqdtingl/Dropbox/Nitz lab Shared/TBG_paper/population_rate_lfp_data_f.mat','files')
load('/home/david/Dropbox/Documents/pubs/bf_theta_beta_gammas/file_all')
load('/home/david/Dropbox/Documents/pubs/bf_theta_beta_gammas/S_markers_times_all')
load('/home/david/Dropbox/Documents/pubs/bf_theta_beta_gammas/F_markers_times_all')
cc=1;
% 
% for r = rstart:8
% 

    list = dir(['/media/david/My Passport/nitzlabdata/LFP_both/'  '*.mat']);
 
    
for i = 1:length(list)
   
         i

   load(['/media/david/My Passport/nitzlabdata/LFP_both/' list(i).name],'ad');
%    ad = add;
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
temp_wave_f=[];
temp_wave_control_f=[];
                
if ~exist(['/media/david/My Passport/nitzlabdata/LFP_ratemaps_both/' list(i).name])  
   if ~isempty(sRec)
   for j = 1:length(ad)
       if ~isnan(ad{j}) & ~isempty(ad{j}) & ad{j} ~= -1 

    
           %% getting amplitude approximate from wavelet
        freq = scal2frq(1:150,'cgau4',1/1000);

        c = cwt(ad{j},freq,'cgau4');
        S = abs(c.*c); clear c;
       wavelet_all = [(100*S./sum(S(:)))];
        phases = angle(hilbert(filtfilt(b,a,ad{j})));
        clear S

                % successful trials
       
                for trial = 1:size(S_markers_times{sRec},1)
                    trial
                    start = ceil((S_markers_times{sRec}(trial,1))*1000);
                    light = ceil((S_markers_times{sRec}(trial,2))*1000);
                    nosepoke = ceil((S_markers_times{sRec}(trial,3))*1000);
                    platecross = ceil((S_markers_times{sRec}(trial,4))*1000);
                    stop = ceil((S_markers_times{sRec}(trial,5))*1000);
                    trialend = ceil((S_markers_times{sRec}(trial,6))*1000);
                    
                    if trialend < size(wavelet_all,2)
                        
                     %% phase/amplitude stuff here
                    pa{j,trial} = zeros(trialend-start,3,35);
%                     for p = start:trialend
%                         for h = 1:35
%                            f = find(round((phases(p-250:p+250)+pi)*100/18)==h);
%                            
%                            pa{j,trial}(p-start+1,1,h) = mean(mean(wavelet_all(20:35,f)'));
%                            pa{j,trial}(p-start+1,2,h) = mean(mean(wavelet_all(45:65,f)'));
%                            pa{j,trial}(p-start+1,3,h) = mean(mean(wavelet_all(80:150,f)'));
%                         end
%                     end
                    for tt = 1%:35
                         add = 0;
                    len =   [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    while length(len) < 240
                    add = add+.01;
                    len =   [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    while length(len) > 240
                    add = add-.01;
                    len =  [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                    add_b = 0;
                    len =   [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    while length(len) < 80
                    add_b = add_b+.01;
                    len =   [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    end
                    while length(len) > 80
                    add_b = add_b-.01;
                    len =  [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    end
                    
                    add_c = 0;
                    len =   [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    while length(len) < 160
                    add_c = add_c+.01;
                    len =   [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    end
                    while length(len) > 160
                    add_c = add_c-.01;
                    len =  [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    end
                    
                    add_d = 0;
                    len =   [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    while length(len) < 60
                    add_d = add_d+.01;
                    len =   [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    end
                    while length(len) > 60
                    add_d = add_d-.01;
                    len =  [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    end
                    
                    
%                     
%                         pa_splined(j,trial,1,:,tt) = ...
%                     [spline(1:size(pa{j,trial}(start+1-start:light+1+1-start),2),pa{j,trial}(start+1-start:light+1+1-start,1,tt),1:size(pa{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2),pa{j,trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa{j,trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2),pa{j,trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(stop+1-start:trialend-start),2),pa{j,trial}(stop+1-start:trialend-start,1,tt),1:size(pa{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{j,trial}(stop+1-start:trialend-start),2))];
%                 
%                     pa_splined(j,trial,2,:,tt) = ...
%                     [spline(1:size(pa{j,trial}(start+1-start:light+1+1-start),2),pa{j,trial}(start+1-start:light+1+1-start,2,tt),1:size(pa{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2),pa{j,trial}(light+1-start:nosepoke+1+1-start,2,tt),1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa{j,trial}(nosepoke+1-start:platecross+1+1-start,2,tt),1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2),pa{j,trial}(platecross+1-start:stop+1+1-start,2,tt),1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(stop+1-start:trialend-start),2),pa{j,trial}(stop+1-start:trialend-start,2,tt),1:size(pa{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{j,trial}(stop+1-start:trialend-start),2))];
%                 
%                 
%                     pa_splined(j,trial,3,:,tt) = ...
%                     [spline(1:size(pa{j,trial}(start+1-start:light+1+1-start),2),pa{j,trial}(start+1-start:light+1+1-start,3,tt),1:size(pa{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2),pa{j,trial}(light+1-start:nosepoke+1+1-start,3,tt),1:size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa{j,trial}(nosepoke+1-start:platecross+1+1-start,3,tt),1:size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2),pa{j,trial}(platecross+1-start:stop+1+1-start,3,tt),1:size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa{j,trial}(stop+1-start:trialend-start),2),pa{j,trial}(stop+1-start:trialend-start,3,tt),1:size(pa{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa{j,trial}(stop+1-start:trialend-start),2))];
%                 
                    end
                    %%
  
                    for tt = 1:150
                    temp_wave(j,trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                    r = randi(5000)-2500;
                    start = ceil((S_markers_times{sRec}(trial,1))*1000)+r;
                    light = ceil((S_markers_times{sRec}(trial,2))*1000)+r;
                    nosepoke = ceil((S_markers_times{sRec}(trial,3))*1000)+r;
                    platecross = ceil((S_markers_times{sRec}(trial,4))*1000)+r;
                    stop = ceil((S_markers_times{sRec}(trial,5))*1000)+r;
                    trialend = ceil((S_markers_times{sRec}(trial,6))*1000)+r;
                    
                    
                      for tt = 1:150
                    temp_wave_control(j,trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
%                     for tt = 1:150
%                     temp_wave2(j,trial,tt,:) = ...
%                     [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:119.5),...
%                     spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:79.5),...
%                     spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:159.5),...
%                     spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:59.5),...
%                     spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:119.5)];
%                     end
                    end

                end
        
                % failure trials
              
               
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
                    pa_f{j,trial} = zeros(trialend-start,3,35);
%                     for p = start:trialend
%                         f=[];
%                         for h = 1:35
%                            f = find(round((phases(p-250:p+250)+pi)*100/18)==h);
%                            
%                            pa_f{j,trial}(p-start+1,1,h) = mean(mean(wavelet_all(20:35,f)'));
%                            pa_f{j,trial}(p-start+1,2,h) = mean(mean(wavelet_all(45:65,f)'));
%                            pa_f{j,trial}(p-start+1,3,h) = mean(mean(wavelet_all(80:150,f)'));
%                            if isempty(f)
%                                pa_f{j,trial}(p-start+1,1,h) = nan;
%                                pa_f{j,trial}(p-start+1,2,h) = nan;
%                                pa_f{j,trial}(p-start+1,3,h) = nan;
%                            end
%                         end
%                     end
                    for tt = 1%:35
                    add = 0;
                    len =   [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    while length(len) < 240
                    add = add+.01;
                    len =   [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    while length(len) > 240
                    add = add-.01;
                    len =  [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                    add_b = 0;
                    len =   [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    while length(len) < 80
                    add_b = add_b+.01;
                    len =   [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    end
                    while length(len) > 80
                    add_b = add_b-.01;
                    len =  [spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2))];
                    end
                    
                    add_c = 0;
                    len =   [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    while length(len) < 160
                    add_c = add_c+.01;
                    len =   [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    end
                    while length(len) > 160
                    add_c = add_c-.01;
                    len =  [spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2))];
                    end
                    
                    add_d = 0;
                    len =   [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    while length(len) < 60
                    add_d = add_d+.01;
                    len =   [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    end
                    while length(len) > 60
                    add_d = add_d-.01;
                    len =  [spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2))];
                    end


                    
                    
%                         pa_splined_f(j,trial,1,:,tt) = ...
%                     [spline(1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2),pa_f{j,trial}(start+1-start:light+1+1-start,1,tt),1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2),pa_f{j,trial}(light+1-start:nosepoke+1+1-start,1,tt),1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start,1,tt),1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2),pa_f{j,trial}(platecross+1-start:stop+1+1-start,1,tt),1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(stop+1-start:trialend-start),2),pa_f{j,trial}(stop+1-start:trialend-start,1,tt),1:size(pa_f{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{j,trial}(stop+1-start:trialend-start),2))];
%                 
%                     pa_splined_f(j,trial,2,:,tt) = ...
%                     [spline(1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2),pa_f{j,trial}(start+1-start:light+1+1-start,2,tt),1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2),pa_f{j,trial}(light+1-start:nosepoke+1+1-start,2,tt),1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start,2,tt),1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2),pa_f{j,trial}(platecross+1-start:stop+1+1-start,2,tt),1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(stop+1-start:trialend-start),2),pa_f{j,trial}(stop+1-start:trialend-start,2,tt),1:size(pa_f{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{j,trial}(stop+1-start:trialend-start),2))];
%                 
%                 
%                     pa_splined_f(j,trial,3,:,tt) = ...
%                     [spline(1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2),pa_f{j,trial}(start+1-start:light+1+1-start,3,tt),1:size(pa_f{j,trial}(start+1-start:light+1+1-start),2)./(119.5+add):size(pa_f{j,trial}(start+1-start:light+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2),pa_f{j,trial}(light+1-start:nosepoke+1+1-start,3,tt),1:size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)./(79.5+add):size(pa_f{j,trial}(light+1-start:nosepoke+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2),pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start,3,tt),1:size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)./(159.5+add):size(pa_f{j,trial}(nosepoke+1-start:platecross+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2),pa_f{j,trial}(platecross+1-start:stop+1+1-start,3,tt),1:size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)./(59.5+add):size(pa_f{j,trial}(platecross+1-start:stop+1+1-start),2)),...
%                     spline(1:size(pa_f{j,trial}(stop+1-start:trialend-start),2),pa_f{j,trial}(stop+1-start:trialend-start,3,tt),1:size(pa_f{j,trial}(stop+1-start:trialend-start),2)./(119.5+add):size(pa_f{j,trial}(stop+1-start:trialend-start),2))];
                
                    end
                    %%
  
                    for tt = 1:150
                    temp_wave_f(j,trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
                    r = randi(10000)-5000;
                    start = ceil((F_markers_times{sRec}(trial,1))*1000)+r;
                    light = ceil((F_markers_times{sRec}(trial,2))*1000)+r;
                    nosepoke = ceil((F_markers_times{sRec}(trial,3))*1000)+r;
                    platecross = ceil((F_markers_times{sRec}(trial,4))*1000)+r;
                    stop = ceil((F_markers_times{sRec}(trial,5))*1000)+r;
                    trialend = ceil((F_markers_times{sRec}(trial,6))*1000)+r;
                    
                        for tt = 1:150
                    temp_wave_control_f(j,trial,tt,:) = ...
                    [spline(1:size(wavelet_all(:,start:light+1),2),wavelet_all(tt,start:light+1),1:size(wavelet_all(:,start:light+1),2)./(119.5+add):size(wavelet_all(:,start:light+1),2)),...
                    spline(1:size(wavelet_all(:,light:nosepoke+1),2),wavelet_all(tt,light:nosepoke+1),1:size(wavelet_all(:,light:nosepoke+1),2)./(79.5+add_b):size(wavelet_all(:,light:nosepoke+1),2)),...
                    spline(1:size(wavelet_all(:,nosepoke:platecross+1),2),wavelet_all(tt,nosepoke:platecross+1),1:size(wavelet_all(:,nosepoke:platecross+1),2)./(159.5+add_c):size(wavelet_all(:,nosepoke:platecross+1),2)),...
                    spline(1:size(wavelet_all(:,platecross:stop+1),2),wavelet_all(tt,platecross:stop+1),1:size(wavelet_all(:,platecross:stop+1),2)./(59.5+add_d):size(wavelet_all(:,platecross:stop+1),2)),...
                    spline(1:size(wavelet_all(:,stop:trialend+1),2),wavelet_all(tt,stop:trialend+1),1:size(wavelet_all(:,stop:trialend+1),2)./(119.5+add):size(wavelet_all(:,stop:trialend+1),2))];
                    end
                    
%                     for tt = 1:150
%                     temp_wave2(j,trial,tt,:) = ...
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
       
   end
%    load(['D:\nitzlabdata/LFP_data_phase_amp/' list(i).name]);
   save(['/media/david/My Passport/nitzlabdata/LFP_ratemaps_both/' list(i).name],'wavelet_all','temp*','pa*','-v7.3')
%   
   
    end
    
% matchedRecs(cc,:) = [i,sRec];
% lfp_rec_wave{cc} = (temp_wave);
% lfp_rec_wave_control{cc} = (temp_wave_control);
% lfp_rec_wave_control_f{cc} = (temp_wave_control_f);
% lfp_rec_pa{cc} = (pa_splined);
% 
% 
% lfp_rec_wave_f{cc} = (temp_wave_f);
% lfp_rec_pa_f{cc} = (pa_splined_f);

clear pa*

cc=1+cc
 end
% end



end


%    save(['/media/nitzlabdata/My Passport/nitzlabdata/LFP_wavelet_ratemaps_BF_PPC.mat'],'-v7.3')

end