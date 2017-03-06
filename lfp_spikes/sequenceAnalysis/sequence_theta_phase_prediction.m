function [aic cellseq logl logl_ind dev] = sequence_theta_phase_prediction(tts,spikes,cellseq,phases,ff)
% cellseq = cell2mat(cellseq);
figure(3)
sp = zeros(length(spikes),length(phases));
for n = 1:length(spikes)
sp(n,round(1000*spikes{n}))=1;
sp(n,:) = smoothts(sp(n,:),'b',50);
end
%% run single spike train GLM's here
for cell = 1:size(sp,1)
   [r_ind dev_ind(cell)] = glmfit(sp(cell,:),phases,'normal','link','identity'); 
   yfit_ind(cell,:) = glmval(r_ind,sp(cell,:),'identity');
   logl_ind(cell) = -nansum(log(normpdf(phases,yfit_ind(cell,:)')));
end



for i=1:size(cellseq,1)
    c = (cellseq(i,:));
    if length(unique(c)) == ff
        v =1;
        sp_seq = zeros(ff,length(phases));
        for j=c
            sp_seq(v,round(spikes{j}*1000)) = 1;
            sp_seq(v,:) = (sp_seq(v,:));
            v=1+v;
        end
%         [ts] = findchain(spikes, c, .02, 100);
        ts = tts{i};
        seq = zeros(1,length(phases));
        seq(round(ts*1000))=1;
       
        [r_spikes dev(i,1)] = glmfit(sp_seq',phases,'normal','link','identity');
        [r_seq dev(i,2)] = glmfit(seq',phases,'normal','link','identity');
        [r_sum dev(i,3)] = glmfit(sum(sp_seq)',phases,'normal','link','identity');

        yfit_spike = glmval(r_spikes,sp_seq','identity');
        yfit_seq = glmval(r_seq,seq,'identity');
        yfit_sum = glmval(r_sum,sum(sp_seq)','identity');
%         
        logl(i,1) = nansum(log(normpdf(phases,yfit_spike)));
        logl(i,2) = nansum(log(normpdf(phases,yfit_seq)));
        logl(i,3) = nansum(log(normpdf(phases,yfit_sum)));
%                 i
        hold on
        aic(i,:) = aicbic([logl(i,:) logl_ind],[ff 1 1 repmat(1,size(sp,1),1)']);
        plot(i*ones(size(sp,1),1),aic(i,4:end)-mean(aic(i,:)),'.b')
        plot([i],aic(i,1)-mean(aic(i,:)),'.m')
        plot([i],aic(i,2)-mean(aic(i,:)),'.g')
        plot([i],aic(i,3)-mean(aic(i,:)),'.r')
        title(['sequence length: ' num2str(ff)])
        
        pause(.1)
    end
end



end


