% used with DT2 27th initially
% 
% 
load('DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226.behavior.mat')
load('DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226.firingMaps.cellinfo.mat','firingMaps')
sessionInfo = bz_getSessionInfo;
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
spikes = bz_GetSpikes;
[b a] = butter(4,[4/625 12/625],'bandpass');
filt = filtfilt(b,a,double(lfp.data));

nCells = length(spikes.times);
nTrials = size(firingMaps.rateMaps{1},2);
% set trial condition here
trials = find(behavior.events.trialConditions==1);

for i=1:nTrials
st = behavior.events.trialIntervals(trials(i),1);
sto = behavior.events.trialIntervals(trials(i),2);
f{i} = filt(round(1250*st):round(1250*sto));
end

% 
% 
ang = angle(hilbert(filt));

 for i=1:nTrials
st = behavior.events.trialIntervals(trials(i),1);
sto = behavior.events.trialIntervals(trials(i),2);
for j=1:length(spikes.times)
sp{i}{j}=spikes.times{j}(intersect(find(spikes.times{j}>st),find(spikes.times{j}<sto)))-st;
end
end
for i=1:nTrials
[a b]=findpeaks(f{i}); b=b./1250;
for j=1:length(b)-1
for k=1:nCells
cycles{i}{j}{k} = sp{i}{k}(intersect(find(sp{i}{k}>b(j)),find(sp{i}{k}<b(j+1))));
angles{i}{j}{k} = ang(ceil(sp{i}{k}(intersect(find(sp{i}{k}>b(j)),find(sp{i}{k}<b(j+1))))*1250));
end
end
end
% 
% % make control cycles here
for j=1:nTrials
    [a b_control{j}]=findpeaks(f{j}); b_control{j}=b_control{j}./length(f{j});
end
% 
for i=1:nTrials
[a b]=findpeaks(-f{i}); b=b./1250;
bb = b./length(f{i}).*1250;
for j=1:length(b)
    cycles_control{i}{j} = cell(1,nCells);
    phase_angles{i}{j} = cell(1,nCells);
%     for k=1:nCells
%         cycles{i}{j}{k} = sp{i}{k}(intersect(find(sp{i}{k}>b(j)),find(sp{i}{k}<b(j+1))));
%     end
    for jj=1:nTrials
        if jj ~= j
        [a near] = min(abs(bb(j)-b_control{jj})); 
        if near == length(b_control{jj})
            near = near -1; % sets back one if at end of trial
        end
        for k=1:nCells
            spk = [intersect(find(sp{jj}{k}>b_control{jj}(near)*length(f{jj})./1250),find(sp{jj}{k}<b_control{jj}(near+1)*length(f{jj})./1250))];
            if ~isempty(spk)
        cycles_control{i}{j}{k} = [cycles_control{i}{j}{k},...
            sp{jj}{k}(spk)'];
        ph = angle(hilbert(f{jj}));
        phase_angles{i}{j}{k} = [phase_angles{i}{j}{k}, ph(ceil(1250*sp{jj}{k}(spk)))'];
            end
            end
        end
    end
end
end
% 
for t=1:nTrials
for k=1:length(cycles{t})
    num = length(cell2mat(cycles{t}{k}'));
    r = randperm(length(cell2mat(cycles_control{t}{k})));
    cyc = [];
    ph=[];
    for i=1:nCells
       cyc = [cyc;cycles_control{t}{k}{i}' repmat(i,length(cycles_control{t}{k}{i}),1)]; 
       ph = [ph;phase_angles{t}{k}{i}'];
    end
    cyc = cyc(r(1:num),:); % matching number of random spikes.times from nearby theta cycles (other trials)
    ph = ph(r(1:num));
    for i=1:nCells
        cyc_new{i} = cyc(find(cyc(:,2)==i),1);
        phh{i}=ph(find(cyc(:,2)==i));
    end
if ~isempty(cyc_new)
[c idx] = sort_firstSpike(cyc_new);
for i=1:nCells
    temp{i}=phh{idx(i)};
end
phh = temp; clear temp
for i=1:nCells
if idx(i)==131
subplot(1,2,1)
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
hold on
subplot(1,2,2)
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'.r')
hold on
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'.r')
else
% plot(c{i},i*ones(length(c{i}),1),'.k')
end

hold on
end;end
end
pause(.01)
end



% 
% for i=1:28
% [a b]=findpeaks(ff{i}); b=b./1250;
% for j=1:length(b)-1
% for k=1:nCells
% cycles2{i}{j}{k} = sp{i}{k}(intersect(find(sp{i}{k}>b(j)),find(sp{i}{k}<b(j+1))));
% end
% end
% end
% 
for t=10:nTrials
for k=1:length(cycles{t})
if ~isempty(cycles{t}{k})
[c idx] = sort_meanSpike(cycles{t}{k});
subplot(2,1,1)
plot(linspace(0,length(f{t})/1250,length(f{t})),f{t})
subplot(2,1,2)
sp=[];
for i=1:nCells
    if ~isempty(c{i})
    sp = [sp;c{i} repmat(idx(i),length(c{i}),1)];
    end
end
for i=1:size(sp,1)
    s{i} = sp(i,1);
end
[c idx] = sort_firstSpike(s); clear s
for i=1:length(idx)
if sp(idx(i),2) == 159
    plot(c{i},i*ones(length(c{i}),1),'.g')
    plot(c{i},i*ones(length(c{i}),1),'og')
elseif sp(idx(i),2) == 63
% plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
    plot(c{i},i*ones(length(c{i}),1),'.m')
    plot(c{i},i*ones(length(c{i}),1),'om')
elseif sp(idx(i),2) == 131
% plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
    plot(c{i},i*ones(length(c{i}),1),'.r')
    plot(c{i},i*ones(length(c{i}),1),'or')
else
    plot(c{i},i*ones(length(c{i}),1),'.k')
end
hold on
end
end
end
pause
clf
end
% % 
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
figure
for t=1:nTrials
for k=1:length(cycles{t})
    num = length(cell2mat(cycles{t}{k}'));
    r = randperm(length(cell2mat(cycles_control{t}{k})));
    cyc = [];
    ph=[];
    for i=1:nCells
       cyc = [cyc;cycles_control{t}{k}{i}' repmat(i,length(cycles_control{t}{k}{i}),1)]; 
       ph = [ph;phase_angles{t}{k}{i}'];
    end
    cyc = cyc(r(1:num),:); % matching number of random spikes.times from nearby theta cycles (other trials)
    ph = ph(r(1:num));
    for i=1:nCells
        cyc_new{i} = cyc(find(cyc(:,2)==i),1);
        phh{i}=ph(find(cyc(:,2)==i));
    end
if ~isempty(cyc_new)
[c idx] = sort_meanSpike(cyc_new);
for i=1:nCells
    temp{i}=phh{idx(i)};
end
phh = temp; clear temp
for i=1:nCells
subplot(1,2,1)
if idx(i)==2
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
hold on
subplot(1,2,2)
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'.r')
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'.r')
hold on
end
end
end
end
pause(.01)
end
pause

%% phew... 

% clf(); 
figure
for t=1:nTrials
for k=1:length(cycles{t})
    num = length(cell2mat(cycles{t}{k}'));
    r = randperm(length(cell2mat(cycles_control{t}{k})));
    cyc = [];
    ph=[];
    for i=1:nCells
       cyc = [cyc;cycles_control{t}{k}{i}' repmat(i,length(cycles_control{t}{k}{i}),1)]; 
       ph = [ph;phase_angles{t}{k}{i}'];
    end
    cyc = cyc(r(1:num),:); % matching number of random spikes.times from nearby theta cycles (other trials)
    ph = ph(r(1:num));
    for i=1:nCells
        cyc_new{i} = cyc(find(cyc(:,2)==i),1);
        phh{i}=ph(find(cyc(:,2)==i));
    end
if ~isempty(cyc_new)
[c idx] = sort_meanSpike(cyc_new);
for i=1:nCells
    temp{i}=phh{idx(i)};
end
phh = temp; clear temp
for i=1:nCells
subplot(2,2,1); hold on
if idx(i) ==63
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
elseif idx(i) == 131
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.g'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'og')
elseif idx(i) == 159
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'+m'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'+m')
elseif idx(i) == 96
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.k'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'ok')
end
subplot(2,2,2); hold on;
if idx(i) ==63
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'.r')
hold on
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'.r')
elseif idx(i) == 131
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'.g')
hold on
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'.g')
elseif idx(i) == 159
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'+m')
hold on
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'+m')
elseif idx(i) == 96
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i},'.k')
hold on
plot(k*ones(length(c{i}),1)./length(cycles{t}),phh{i}+2*pi,'.k')
end
end
[c idx] = sort_meanSpike(cycles{t}{k});
for i=1:nCells
subplot(2,2,3); hold on;
if idx(i) ==63
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.r'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'or')
elseif idx(i) == 131
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.g'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'og')
elseif idx(i) == 159
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'+m'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'+m')
elseif idx(i) == 96
plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'.k'); plot(k*ones(length(c{i}),1)./length(cycles{t}),i*ones(length(c{i}),1)./sum(~cellfun(@isempty,c)),'ok')
end
subplot(2,2,4); hold on
ph = angle(hilbert(f{t})); 

if idx(i) ==63
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250)),'.r');
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250))+2*pi,'.r');
elseif idx(i) == 131
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250)),'.g');
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250))+2*pi,'.g');
elseif idx(i) == 159
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250)),'+m');
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250))+2*pi,'+m');
elseif idx(i) == 96
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250)),'.k');
plot(k*ones(length(c{i}),1)./length(cycles{t}),ph(ceil(c{i}*1250))+2*pi,'.k');
end
end
end
end
pause(.01)
t
end

subplot(2,2,1)
title('cycle shuffled theta sequences')
subplot(2,2,2)
title('cycle shuffled precession')
subplot(2,2,3)
title('actual theta sequences')
subplot(2,2,4)
title('actual phase precession')



