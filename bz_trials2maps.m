function [map mapping] = bz_trials2maps(trials)

c=1;
for tt = 1:length(trials)    
    map{tt}=[];
    t_conc=[];
    for t = 1:length(trials{tt})
        t_conc = [trials{tt}{t}(:,:),20*(trials{tt}{t}(:,1)-trials{tt}{t}(1,1))];
    if length(t_conc)>=200
        while length(t_conc)>200+1
        di = pdist(t_conc);
        s = squareform(di);
        s(find(eye(size(s))))=nan;
        [a b] = min(s(:));
        [coords blah] = find(s==a);
        t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
        t_conc(coords(2),:) = [];
        % debug
%         scatter(t_conc(:,1),t_conc(:,2))
%         pause(.01)
        end
        if size(t_conc,1) < 201
           for  j=1:size(t_conc,2)
               t_conc_temp(:,j) = makeLength(t_conc(:,j),201); 
           end
           t_conc = t_conc_temp; clear t_conc_temp
        end
    t_conc_all(t,:,:) = t_conc;
    else
        while length(t_conc)>201
        di = pdist(t_conc);
        s = squareform(di);
        s(find(eye(size(s))))=nan;
        [a b] = min(s(:));
        [coords blah] = find(s==a);
        t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
        t_conc(coords(2),:) = [];
        % debug
        %         scatter(t_conc(:,1),t_conc(:,2))
        %         pause(.01)
        end
        if size(t_conc,1) < 201
        for  j=1:size(t_conc,2)
        t_conc_temp(:,j) = makeLength(t_conc(:,j),201); 
        end
        t_conc = t_conc_temp; clear t_conc_temp
        end
        t_conc_all(t,:,:) = t_conc;
        
    end
    end
    
    if length(trials{tt})>0   
        map{tt} = squeeze(median(t_conc_all(:,:,:),1));
    end


    clear t_conc_all

    disp('finding mapping...')
%     map{tt} = [w,ww];
    for t =1:length(trials{tt})  % all trial types (rotations)
        for p = 1:length(trials{tt}{t})
            [a b] = min(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
                trials{tt}{t}(p,8)-map{tt}(:,8),...
                trials{tt}{t}(p,9)-map{tt}(:,9),...
                trials{tt}{t}(p,10)-map{tt}(:,10),...
                (trials{tt}{t}(p,1)-trials{tt}{t}(1,1))*50-map{tt}(:,1),...  % penalty for time differences
                40*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'));     % penalty for order differences
            mapping{tt}{t}(p,:) = [map{tt}(b,1:end) b trials{tt}{t}(p,1)];
%             plot(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
%                 trials{tt}{t}(p,2)-map{tt}(:,2),...
%                 trials{tt}{t}(p,3)-map{tt}(:,3),...
%                 trials{tt}{t}(p,4)-map{tt}(:,4),...
%                 (trials{tt}{t}(p,5)-trials{tt}{t}(1,5))*20-map{tt}(:,5),...  % penalty for time differences
%                 50*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'))
%             pause
        end
    end
end
end
