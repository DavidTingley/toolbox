function [aic cellseq logl logl_ind dev] = sequence_position_prediction(tts,spikes,cellseq,ff, position)
% cellseq = cell2mat(cellseq);
figure(2)
% if size(position,2) > 2
%    position  = position(:,4,5); 
% end

position(position < 20) = nan;
position = position(:,2);%-position(:,3);
% position(position < -50) = nan;
% position(position > 420) = nan;
sp = zeros(length(spikes),length(position));
for n = 1:length(spikes)
sp(n,round(1000*spikes{n}))=1;
sp(n,:) = smoothts(sp(n,:),'b',50);
end

%% interpolate position to neural sampling frequency here
x = 1:length(position);
xx = 0:(length(position)./length(sp)):length(position)-.01;
p(:,1) = interp1(x,position,xx);
% p(:,2) = interp1(x,position(:,5),xx);
position = p;
%


%% run single spike train GLM's here
for cell = 1:size(sp,1)
   [r_ind dev_ind(cell,1)] = glmfit(sp(cell,:),position(:,1),'normal','link','identity'); 
   yfit_ind(cell,:) = glmval(r_ind,sp(cell,:),'identity');
   n =normpdf(position(:,1),yfit_ind(cell,:)');
   n(n==0) = nan;
   logl_ind(cell,1) = -nansum(log(n));
%    [r_ind dev_ind(cell,2)] = glmfit(sp(cell,:),position(:,2),'normal','link','identity'); 
%    yfit_ind(cell,:) = glmval(r_ind,sp(cell,:),'identity');
%    n =normpdf(position(:,2),yfit_ind(cell,:)');
%    n(n==0) = nan;
%    logl_ind(cell,2) = -nansum(log(n));
end



for i=1:size(cellseq,1)
    c = (cellseq(i,:));
    if length(unique(c)) == ff
        v =1;
        sp_seq = zeros(ff,length(position));
        for j=c
            sp_seq(v,round(spikes{j}*1000)) = 1;
            sp_seq(v,:) = (sp_seq(v,:));
            v=1+v;
        end
%         [ts] = findchain(spikes, c, .02, 100);
%         if size(ts,1) > 1 & iscell(ts(1))
        ts = tts{i};
%         end
        seq = zeros(1,length(position));
        seq(round(ts*1000))=1;
       
        [r_spikes dev(i,1)] = glmfit(sp_seq',position(:,1),'normal','link','identity');
        [r_seq dev(i,2)] = glmfit(seq',position(:,1),'normal','link','identity');
        [r_sum dev(i,3)] = glmfit(sum(sp_seq)',position(:,1),'normal','link','identity');
        
        yfit_spike(1,:) = glmval(r_spikes,sp_seq','identity');
        yfit_seq(1,:) = glmval(r_seq,seq,'identity');
        yfit_sum(1,:) = glmval(r_sum,sum(sp_seq)','identity');
        
        n_spike =normpdf(position(:,1),yfit_spike(1,:)');
        n_spike(n_spike==0) = nan;
        n_seq =normpdf(position(:,1),yfit_seq(1,:)');
        n_seq(n_seq==0) = nan;
        n_sum =normpdf(position(:,1),yfit_sum(1,:)');
        n_sum(n_sum==0) = nan;
        
        logl(i,1,1) = -nansum(log(n_spike));
        logl(i,1,2) = -nansum(log(n_seq));
        logl(i,1,3) = -nansum(log(n_sum));
        
%    _theta
%         [r_spikes dev_spikes(i,2)] = glmfit(sp_seq',position(:,2),'normal','link','identity');
%         [r_seq dev_seq(i,2)] = glmfit(seq',position(:,2),'normal','link','identity');
%         [r_sum dev_sum_spikes(i,2)] = glmfit(sum(sp_seq)',position(:,2),'normal','link','identity');
% 
%         yfit_spike(2,:) = glmval(r_spikes,sp_seq','identity');
%         yfit_seq(2,:) = glmval(r_seq,seq,'identity');
%         yfit_sum(2,:) = glmval(r_sum,sum(sp_seq)','identity');
%         
%         n_spike =normpdf(position(:,2),yfit_spike(2,:)');
%         n_spike(n_spike==0) = nan;
%         n_seq =normpdf(position(:,2),yfit_seq(2,:)');
%         n_seq(n_seq==0) = nan;
%         n_sum =normpdf(position(:,2),yfit_sum(2,:)');
%         n_sum(n_sum==0) = nan;
%         _theta
%         logl(i,2,1) = -nansum(log(n_spike));
%         logl(i,2,2) = -nansum(log(n_seq));
%         logl(i,2,3) = -nansum(log(n_sum));
% %                 i
        
        aic(i,1,:) = aicbic([squeeze(logl(i,1,:)); logl_ind(:,1)]',[ff 1 1 repmat(1,size(sp,1),1)']);
%         aic(i,2,:) = aicbic([squeeze(logl(i,2,:)); logl_ind(:,2)]',[ff 1 1 repmat(1,size(sp,1),1)']);
        
        subplot(2,1,1)
        hold on
        plot(i*ones(size(sp,1),1),squeeze(aic(i,1,4:end))-mean(aic(i,1,:)),'.b')
        plot([i],squeeze(aic(i,1,1)./mean(aic(i,1,:))),'.m')
        plot([i],squeeze(aic(i,1,2)./mean(aic(i,1,:))),'.g')
        plot([i],squeeze(aic(i,1,3)./mean(aic(i,1,:))),'.r')
        subplot(2,1,2)
        hold on
%         plot(i*ones(size(sp,1),1),squeeze(aic(i,1,4:end))-mean(aic(i,1,:)),'.b')
        plot([i],squeeze(aic(i,1,1))./abs(mean(aic(i,1,1))),'.m')
        plot([i],squeeze(aic(i,1,2))./abs(mean(aic(i,1,1))),'.g') %./abs((aic(i,1,1)))
        plot([i],squeeze(aic(i,1,3))./abs(mean(aic(i,1,1))),'.k')
        title(['x position, sequence(g), trains(m), sum trains(k) length: ' num2str(ff)])
        
%         subplot(2,1,2)
%         hold on
%         plot(i*ones(size(sp,1),1),squeeze(aic(i,2,4:end)-mean(aic(i,2,:))),'.b')
%         plot([i],squeeze(aic(i,2,1)-mean(aic(i,2,:))),'.m')
%         plot([i],squeeze(aic(i,2,2)-mean(aic(i,2,:))),'.g')
%         plot([i],squeeze(aic(i,2,3)-mean(aic(i,2,:))),'.r')
%         title(['y position, sequence length: ' num2str(ff)])
        
        pause(.1)
    end
end



end


