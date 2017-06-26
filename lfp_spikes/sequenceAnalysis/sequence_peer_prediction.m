function [aic cellseq logl dev] = sequence_peer_prediction(ts,spikes,cellseq,seqLength)
    figure(1)
    maxlen = 0;
    for n = 1:length(spikes)
        if ~isempty(spikes{n})
        if maxlen < spikes{n}(end)
            maxlen = spikes{n}(end);           
        end
        end
    end
    
    sp = zeros(n,round(maxlen*1000)+10);
    for n = 1:length(spikes)
        if ~isempty(spikes{n})
            sp(n,ceil(spikes{n}*1000))=1;
            sp(n,:) = fastrms(sp(n,:),50);
        end
    end
    count=1;
    for i=1:size(cellseq,1)
        c = cellseq(i,:);
        hold on
        if length(unique(c)) == seqLength
            v = 1;
            sp_seq = zeros(seqLength,round(maxlen*1000)+10);
            for j=c
                sp_seq(v,ceil(spikes{j}*1000)) = 1;
                sp_seq(v,:) = fastrms(sp_seq(v,:),50);
                v=1+v;
            end
            seq = zeros(1,round(maxlen*1000)+10);
            seq(round(ts{i}(:,2)*1000))=1;
            seq = fastrms(seq,50);

            sp_cut = sp; sp_cut(c,:) = [];  % array of spike trains not in the sequence
            for ii = 1:length(spikes)-seqLength
                [r_spikes dev(i,ii,1)] = glmfit(sp_seq',sp_cut(ii,:),'normal','link','identity');
                [r_seq dev(i,ii,2)] = glmfit(seq,sp_cut(ii,:),'normal','link','identity');
                [r_sum dev(i,ii,3)] = glmfit(sum(sp_seq)',sp_cut(ii,:),'normal','link','identity');
                
                
                yfit_spike = glmval(r_spikes,sp_seq','identity');
                yfit_seq = glmval(r_seq,seq,'identity');
                yfit_sum = glmval(r_sum,sum(sp_seq)','identity');
                %         
                p = poisspdf(sp_cut(ii,:),yfit_spike');
                p(p==0)=nan;
                logl(i,ii,1) = nansum(log(p));
                p = poisspdf(sp_cut(ii,:),yfit_seq');
                p(p==0)=nan;
                logl(i,ii,2) = nansum(log(p));
                p =  poisspdf(sp_cut(ii,:),yfit_sum');
                p(p==0)=nan;
                logl(i,ii,3) = nansum(log(p));
        
        
                aic(i,ii,:) = aicbic([logl(i,ii,:)],[4 1 1]);
                plot([count],squeeze(aic(i,ii,1))./mean(aic(i,ii,:)),'.k')
                plot([count],squeeze(aic(i,ii,2))./mean(aic(i,ii,:)),'.g')
                plot([count],squeeze(aic(i,ii,3))./mean(aic(i,ii,:)),'.r')
                count=1+count;
            end
%             plot([i],mean(squeeze(aic(i,:,1)))','.k')
%             plot([i],mean(squeeze(aic(i,:,2)))','.g')
%             plot([i],mean(squeeze(aic(i,:,3)))','.r')
            title('AIC scores: green-seq, red-sum, black-all trains')
            i
            pause(.1)
        end
    end
    
end