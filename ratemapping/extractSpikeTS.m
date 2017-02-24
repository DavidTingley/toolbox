function tfiles = extractSpikeTS(plxsfilename,varargin)
% extractSpikeTS - reads a plx file and extracts timestamps for all spikes,
% separated by channel and by unit. 
%   input: a string containing the filename of the plx file to be read.
%      Optional second input: a string containing the full path to a folder
%      in which to save the extracted t-files. If second input not
%      provided, then the files will be saved in the folder containing the
%      given .plx file. 
%   output: a struct containing the tfile_list and all sig-files as
%       separate fields; these are also all saved as separate files in a
%       new folder called 'tfiles', created in the given save-folder (or
%       the plx-containing folder, if only one input provided). 
% note: for some reason, the generated tfile_list seems to display
% improperly in Notepad, lacking all line-returns. However, it displays
% correctly in other txt-file-reading programs, and the ratemapping program
% seems to read it correctly. 


%help plx_info
%{
  plx_info(filename, fullread) -- read and display .plx file info
 
  [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
 
  INPUT:
    filename - if empty string, will use File Open dialog
    fullread - if 0, reads only the file header
               if 1, reads the entire file
 
  OUTPUT:
    tscounts - 2-dimensional array of timestamp counts for each unit
       tscounts(i, j) is the number of timestamps for channel j-1, unit i
                                 (see comment below)
    wfcounts - 2-dimensional array of waveform counts for each unit
      wfcounts(i, j) is the number of waveforms for channel j-1, unit i
                                 (see comment below)
    evcounts - 1-dimensional array of external event counts
      evcounts(i) is the number of events for event channel i
 
    contcounts - 1-dimensional array of sample counts for continuous channels
      contcounts(i) is the number of continuous for slow channel i-1
 
  Note that for tscounts, wfcounts, the unit,channel indices i,j are off by one. 
  That is, for channels, the count for channel n is at index n+1, and for units,
   index 1 is unsorted, 2 = unit a, 3 = unit b, etc
  The dimensions of the tscounts and wfcounts arrays are
    (NChan+1) x (MaxUnits+1)
  where NChan is the number of spike channel headers in the plx file, and
  MaxUnits is 4 if fullread is 0, or 26 if fullread is 1. This is because
  the header of a .plx file can only accomodate 4 units, but doing a
  fullread on the file may show that there are actually up to 26 units
  present in the file. Likewise, NChan will have a maximum of 128 channels
  if fullread is 0.
  The dimension of the evcounts and contcounts arrays is the number of event
  and continuous (slow) channels. 
  The counts for slow channel 0 is at contcounts(1)
%}
%help plx_ts
%{
  plx_ts(filename, channel, unit): Read spike timestamps from a .plx file
 
  [n, ts] = plx_ts(filename, channel, unit)
 
  INPUT:
    filename - if empty string, will use File Open dialog
    channel - 1-based channel number
    unit  - unit number (0- unsorted, 1-4 units a-d)
 
  OUTPUT:
    n - number of timestamps
    ts - array of timestamps (in seconds)
%}

tfiles = []; 

[inpath, infname, inext] = fileparts(plxsfilename); 
if(isempty(inpath))
    % if the input includes only a relatve path, we use the current dir. 
    plxspath = pwd; 
    plxsfilename = fullfile(plxspath, plxsfilename);
else
    % otherwise, use the input path and filename as given. 
    plxspath = inpath; 
end

origdir = pwd; 
% this line creates an object that will execute the specified clean-up
% function when extractSpikeTS stops running, regardless of the reason. 
cObj = onCleanup(@() cd(origdir)); 

% handle user preferences regarding the output save-directory. 
if(~isempty(varargin))
    outsavedir = varargin{1}; 
else
    outsavedir = plxspath; 
end


% extract a list of how many waveform-timestamps were recorded for each
% unit, on each channel. 
spikeTScounts = plx_info(plxsfilename,1);

% process 'spikeTScounts' to eliminate all empty units and empty channels. 
%
% these are the channel numbers that contain some non-zero number of spikes. 
spikechannums = find(any(spikeTScounts,1)); 
% and these are the full channels, containing all units for each. 
spikechans = spikeTScounts(:,spikechannums); 

% we cut down our unit-lists to the length of the largest channel. 
chanunitvec = cell(1,length(spikechannums)); 
numunits = zeros(size(chanunitvec)); 
for ii=1:length(spikechannums)
    chanunitvec{ii} = find(spikechans(:,ii)); 
    numunits(ii) = length(chanunitvec{ii}); 
end


spikechans = spikechans(1:max(numunits),:);
% since we're working with stereotrodes or tetrodes, some channels will be
% identical. Let's find and eliminate all duplicate channels. 
uniquespikechanvec = true(1,size(spikechans,2)); 
for ii=2:size(spikechans,2)
   if(isequal(spikechans(:,ii),spikechans(:,ii-1)))
       uniquespikechanvec(ii) = false; 
   end
end
numunits = numunits(uniquespikechanvec); 
uniquespikechannums = spikechannums(uniquespikechanvec); 
uniquespikechans = spikechans(:,uniquespikechanvec); 


% let's also remove all 'unsorted' units, and exclude any channels that had
% spikes but no other, sorted units.
for ii=1:length(numunits)
    if(numunits(ii)<=1)
        uniquespikechannums(ii) = 0; 
        uniquespikechans(:,ii) = 0; 
    end
end
uniquespikechannums = uniquespikechannums(uniquespikechannums~=0); 
uniquespikechans(1,:) = [];
uniquespikechans = uniquespikechans(:,any(uniquespikechans,1)); 
numunits = numunits-1; 
numunits = numunits(numunits~=0); 

% if, after all that processing, we're left with no channels and no units,
% return an error -- ask the user to supply us with a CUT plx file. 
if(isempty(numunits))
    error('No units found in file: ''%s''. \nEnsure that the given .plx file has been cluster-cut. ',plxsfilename); 
end

% construct unique names for each channel corresponding to the
% 'sig001'-style names they would have if we exported them via NEX. 
sigs = struct('unitcounts',[uniquespikechannums',numunits'],'tscounts',[]); 
signames = cell(length(uniquespikechannums),1); 
for ii=1:length(uniquespikechannums)
    % the channel numbers returned by our plexon utilities are off-by-one
    % from their displayed 'sig' numbers. 
    channumstr = num2str(uniquespikechannums(ii)-1); 
    switch length(channumstr)
        case 1
            signames{ii} = ['sig00', channumstr]; 
        case 2
            signames{ii} = ['sig0', channumstr]; 
        case 3
            signames{ii} = ['sig', channumstr]; 
    end
    
    sigs.(signames{ii}) = []; 
end

% just for fun (and backward-compatibility), let's also save these as
% separate 'tfiles', including a 'tfilelist'.
tfiles = struct('tfile_list',[]); 

% now, let's extract the actual spike timestamps for each of the units
% we've identified. 
sigs.tscounts = cell(size(sigs.unitcounts)); 
for ii=1:length(uniquespikechannums)
    sigs.(signames{ii}) = cell(1,numunits(ii)); 
    sigs.tscounts{ii,1} = signames{ii}; 
    sigs.tscounts{ii,2} = zeros(numunits(ii),1); 
    for ij=1:numunits(ii)
        %disp([signames{ii},' unit ',num2str(ij)])
        [sigs.tscounts{ii,2}(ij),sigs.(signames{ii}){ij}] = plx_ts(plxsfilename,uniquespikechannums(ii),(ij));
        % now the tfiles  (char(97) == 'a')
        tfiles.([signames{ii},char(ij+96)]) = sigs.(signames{ii}){ij}; 
    end
end

tfiles.tfile_list = fieldnames(tfiles); 
% remove 'tfile_list' from the tfile_list field. 
tfiles.tfile_list(1) = []; 


% now we save the results, in either the user-supplied save folder or the
% directory where we found the plx file. Either way, we place our output
% files in a new subfolder called 'tfiles'. 
cd(outsavedir); 
% if(~exist(fullfile(pwd, 'tfiles'),'dir'))
%     mkdir('tfiles');
% end
% cd('tfiles'); 
% save each of the t-files as a separate ASCII file containing a list of
% all spike timestamps for that unit, exactly as if exported from NEX. 
for ii=1:length(tfiles.tfile_list)
    save(tfiles.tfile_list{ii}, '-struct', 'tfiles', tfiles.tfile_list{ii}, '-ascii'); 
end

try
    tflistFID = fopen('tfile_list.txt', 'w+');
    for ii=1:(length(tfiles.tfile_list)-1)
        fprintf(tflistFID, '%s\n', tfiles.tfile_list{ii});
    end
    fprintf(tflistFID, '%s', tfiles.tfile_list{end});
    %save('tfile_list.txt', '-struct', 'tfiles', 'tfile_list', '-ascii');
catch ME
    fclose(tflistFID); 
    error(ME); 
end

fclose(tflistFID); 

end
