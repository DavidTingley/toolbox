%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO DO
% fix binning procedure so they overlap instead of border eachother
% find a valid predictability/performance variable

% this file creates the spike trains used for the GLM model

%% create PETHs with varying/overlapping time-bins for all cells
sigs = '/media/nitzlab/data/Dropbox/Dropbox/Spike_trains/';

% load('C:\Users\Nitz Lab 2\Dropbox_tingley\Dropbox\Dropbox\num_cells.mat')
load('/media/nitzlab/data/Dropbox/Dropbox/num_cells_new')
ncells = num_cells;
rec_start = 40;

if rec_start > 1 
    cell_num_total = sum(ncells(1:rec_start-1));
elseif rec_start == 1
    cell_num_total = 0;
end


for rec = rec_start

wire = [];
bad_rec = 0;
      

for cell_num= 1:ncells(rec)

    cell_num_total = cell_num_total+1
  
        
if ncells(rec) > 1
    if cell_num ==1
    load([sigs 'SIG_' num2str(cell_num_total,'%05d') '.mat'],'spikes_all',...
        'tfile_list','success_markers','recording','ratname')
    end
    if cell_num > 1
    load([sigs 'SIG_' num2str(cell_num_total,'%05d') '.mat'],'spikes_all',...
        'tfile_list')        
    end

     if  length(spikes_all) < success_markers(end,5)*1000+6000

        spikes_all  = [spikes_all; zeros(round(success_markers(end,5)*1000+6000-length(spikes_all)),1)];

     end

    %%%%%%%%%%%%%%%%%%%%  BUILD IN wire seperator

    if cell_num == 1  % only needs to run on the first iteration through recording
        wire = zeros(length(tfile_list),1);
        parfor t= 1:length(tfile_list)
            wire(t) = str2double(tfile_list{t}(5:6));
        end
    end 


    %%%%%%%%%%%%%%%%%%%%
    for i = 1:length(success_markers)
         success_trains(i,:) = test_func2_pairwise(success_markers,spikes_all,i);
    end
    s{cell_num} = success_trains;
    
clear  success_trains spiral* tt s_all*


end 
end 

    if ncells(rec) > 1 & bad_rec ~= 1 & length(wire) == size(s,2)
    disp(['Saving files for ' ratname '_rec_' recording '...'])
    
%     if ~exist(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise/' ratname ...
%         '/' num2str(recording)])  % checks if data has already been generated for recording
    mkdir(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise/' ratname ...
        '/' num2str(recording)])
    
    parfor i = 0:600 % ms window
         [s_all s_all_norm] = paradd_pairwise(ncells,rec,i,s);
          fsave(['/media/nitzlab/Iomega HDD/David/recording_trains_pairwise/' ...
              ratname '/' num2str(recording) '/' ratname '_rec_' recording...
            '_winsize_' num2str(i)],s_all,s_all_norm,wire,recording,ratname,ncells,rec);
    end
    clearvars -except rec sigs cell_num_total ncells temp
    end
%     end

end





