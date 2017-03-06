function [cycles] = findThetaCycles(lfp,varargin)

% Count the number of spikes per cycle in the ongoing oscillatory LFP
% (e.g. during theta).
%
%  USAGE
%
%    [cycles] = findThetaCycles(lfp,samplingRate,varargin)
%
%    
%    lfp         LFP single channel
%
%  OUTPUT
%
%    cycles         list of [start,stop] samples for each cycle
%

% this function takes a single channel with theta oscillations and returns the peak times for each cycle

% check that lfp is a single vector
if min(size(lfp)) > 1
	error('LFP variable doesnt seem to be a single column..')
end

% default args
freqRange = [6 12];
samplingRate = 1250; 

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help findThetaCycles">findThetaCycles</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		% case 'show',
		% 	show = varargin{i+1};
		% 	if ~isstring_FMAT(show,'off','on'),
		% 		error('Incorrect value for property ''show'' (type ''help <a href="matlab:help findThetaCycles">findThetaCycles</a>'' for details).');
		% 	end
		case 'freqrange',
			freqRange = varargin{i+1};
			% if ~isdmatrix(freqRange)
			% 	error(['Parameter ' num2str(i+2) ' is not a property.'])
			% end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help findThetaCycles">findThetaCycles</a>'' for details).']);
	end
end

nyquistLimit = samplingRate/2;

[b a] = butter(4,[freqRange(1)/nyquistLimit freqRange(2)/nyquistLimit],'bandpass');
filtLFP = filtfilt(b,a,lfp);
powerLFP = fastrms(filtLFP,480*samplingRate/1000);
anglLFP = angle(hilbert(filtLFP));

[pks locs]=findpeaks(sin(anglLFP));

cycles = [0 locs'; locs' length(anglLFP)]';
end