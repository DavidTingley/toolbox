 clear 
rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
p = pwd;

                
ratwire_loc_list_parietal = [1:16, zeros(1,length(17:32)),33:48; ...  % AA8
    1:48; ...                                                   % DN3
    1:32,zeros(1,length(33:48));  ...                          % NS1
    1:16,zeros(1,length(17:32)),33:48; ...                    % NS2
    1:16,zeros(1,length(17:32)),33:48; ...                     % NS3
    1:16,zeros(1,length(17:32)),33:48; ...                     % NS4
    1:16,zeros(1,length(17:48));                             % NS5
    1:16,zeros(1,length(17:48));                             % NS8
    1:16,zeros(1,length(17:48))];                             % NS9
                    % this lists the wires that ARE NOT in the Pcx for each rat
                    % order is AA8 DN3 NS1 NS2 NS3 NS4

% getting rats directory
for abc = 1:length(rats)
    ratpath = [p '/' sprintf(rats(abc,:)) '/'];
    
    % getting recordings directory
    recordings = dir([ratpath '/*']);
    recordings = recordings(3:length(recordings));
    
    
    % recording check and make list of folders that need files
    [list_to_ratemap,list_to_ratemap_name, list_need_files] = rate_map_check(ratpath,recordings);

    
    % ratemap recordings with all files 
for rr = 1:length(list_to_ratemap)
        r = list_to_ratemap{rr};
        spiral_info = dir([ratpath recordings(r).name '/' '*spiral*']);
        dvt = dir([ratpath recordings(r).name '/*1.DVT']);
        dvt2 = dir([ratpath recordings(r).name '/*2.DVT']);
        
        spr = distInfo(rats(abc,:),recordings(r).name);
        
        
        connect_issue = dir([ratpath recordings(r).name '/*reverse*']);
        if ~isempty(connect_issue)
            disp('connectors may be reversed')
            disp([ratpath recordings(r).name])
            if 'E:\nitzlabdata/NS3/38' ~= [ratpath recordings(r).name]
            return
            end
            ratwire_loc_list_parietal(5,:) = [1:16,17:32,zeros(1,length(33:48))];
        end
        if isempty(connect_issue)
            % reset ratwire_loc_list variable
            ratwire_loc_list_parietal(5,:) = [1:16,zeros(1,length(17:32)),33:48];
        end 
         
%         % Head direction
                d = load([ratpath recordings(r).name '/' dvt.name]);
        if size(d,2) == 6
              [headDir] = dvtHeadDir([ratpath recordings(r).name '/' dvt.name]);
%               [headDir] = getHeadDirs([ratpath recordings(r).name '/' dvt.name]); 
        end
        if size(d,2) < 6
            d = load([ratpath recordings(r).name '/' dvt2.name]);
            [headDir] = dvtHeadDir([ratpath recordings(r).name '/' dvt2.name]);
%             [headDir] = getHeadDirs([ratpath recordings(r).name '/' dvt2.name]); 
        end    
             
         
        
        % gets actual vs went to files and keystrokes
%           if r == list_to_ratemap{rr}
%         [status,avwt] = fetchAVWT(rats(abc).name,recordings(r).name,[ratpath recordings(r).name '/' dvt.name]);
%           end
if abc < 9
        load([ratpath recordings(r).name '/' rats(abc,:) '_rec' recordings(r).name...
            '_spiral_info.mat'])
elseif abc ==9
        load([ratpath recordings(r).name '/' rats(abc,:) '_rec' recordings(r).name...
            '_NP_spiral_info.mat'])
end
        load([ratpath recordings(r).name '/' rats(abc,:) '_rec' recordings(r).name...
            '_platecross_and_stoptimes.mat'])
        pos = d(:,2:6); 
        [tfile_list]=textread([ratpath recordings(r).name '/tfile_list.txt'],'%s');
        
        %%%% have to pull out parietal wires here
        
        % the following code finds sig files that match parietal wires from
        % the rat_wire_loc_list variable using the tfile list and removes
        % them from the tfile_list variable, so that they are not
        % ratemapped.  Eventually they will be ratemapped seperately for
        % parietal analysis
        count = 1;
        for i = 1:length(tfile_list)
           for j = 1:length(ratwire_loc_list_parietal(abc,:))
               if ratwire_loc_list_parietal(abc,j) ~= 0
                 wire = findstr(tfile_list{i},num2str(ratwire_loc_list_parietal(abc,j),'%03d'));
                 if wire ~= 0 
                    deletions(count) = i;
                    count = 1 + count;
                 end
               end
           end

        end   
        if count ~= 1  
        tfile_list(deletions) = [];
        end
        clear deletions
        
        %%% END tfile_list editing  for only BF cells
        if isempty(tfile_list) == 0  % if only cells were parietal, skip ratemapping
        tfile_list_char(1:length(tfile_list),:)=char(tfile_list);
        avwt=[];
        datapath = [ratpath recordings(r).name '/'];
        t = dir([datapath '*.mat']);
        for i = 1:length(t)
            if ~isempty(regexpi(t(i).name,'.*avwt.*'))
                avwt = load([datapath t(i).name]);
            end
        end
        %{
        rats(abc).name = 'AA8';
        recordings(r).name = '04'; 
        ratpath = [p '/' sprintf(rats(abc).name) '/'];
        datapath = [ratpath recordings(r).name '/'];
        %}
 
%         hd = processArenaData(rats(abc).name,recordings(r).name,datapath);
        hd = [];
            t = dir([datapath '/*.mat']);
            for i = 1:length(t)
               n = regexpi(t(i).name,'.*(light|final).*','match','once'); 
               if ~isempty(n)
                   load([datapath n])
                   hd = 1;
                   break
               end
            end
            disp(datapath)
        if ~isempty(hd)
            % should be markers 
            list = who; 
            for i = 1:length(list)
            if length(list{i}) == length('success_beg')
            if list{i} == 'success_beg'
            s = eval(list{i});
            end
            end
            if length(list{i}) == length('success_begin')
            if list{i} == 'success_begin'
            s = eval(list{i});
            end
            end
            if length(list{i}) == length('correct_begin')
            if list{i} == 'correct_begin'
            s = eval(list{i});
            end
            end
            end
            for i = 1:length(list)
            if length(list{i}) == length('fail_beg')
            if list{i} == 'fail_beg'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('fail_begin')
            if list{i} == 'fail_begin'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('incorrect_begin')
            if list{i} == 'incorrect_begin'
            f = eval(list{i});
            end
            end
            end
            for i = 1:length(list)
            if length(list{i}) == length('success_end')
            if list{i} == 'success_end'
            s_end = eval(list{i});
            end
            end
            if length(list{i}) == length('success_end')
            if list{i} == 'success_end'
            s_end = eval(list{i});
            end
            end
            if length(list{i}) == length('correct_end')
            if list{i} == 'correct_end'
            s_end = eval(list{i});
            end
            end
            end
            for i = 1:length(list)
            if length(list{i}) == length('fail_end')
            if list{i} == 'fail_end'
            f_end = eval(list{i});
            end
            end
            if length(list{i}) == length('fail_end')
            if list{i} == 'fail_end'
            f_end = eval(list{i});
            end
            end
            if length(list{i}) == length('incorrect_end')
            if list{i} == 'incorrect_end'
            f_end = eval(list{i});
            end
            end
            end
            
            clear *arker* *end *beg* *DVT* *anogo*
              if ~isempty(avwt)
            n = num2str(find(spiral_events(:,3)==7)/3);
            a = fieldnames(avwt);
            for temp = 1:length(n)
                avwt_succ(temp,:) = eval(['avwt.' a{:} '(' num2str(n(temp,:)) ',:)']);
            end
            n = num2str(find(spiral_events(:,3)==4)/3);
            for temp = 1:length(n)
                avwt_fail(temp,:) = eval(['avwt.' a{:} '(' num2str(n(temp,:)) ',:)']);
            end
              
   
            
            
            
            if size(s,1) == length(spiral_events(find(spiral_events(:,3)==9),2))...
                    && size(f,1) == length(spiral_events(find(spiral_events(:,3)==10),2))
                disp(datapath)

            success_markers_times(:,2) = s(:,4);  % light flashes
            failure_markers_times(:,2) = f(:,4);
            
            success_markers_times(:,3) = spiral_events(find(spiral_events(:,3)==9),2); % nose-pokes
            failure_markers_times(:,3) = spiral_events(find(spiral_events(:,3)==10),2);
            
            if abc == 9
                success_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==9),2)/60; % plate crosses
                failure_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==10),2)/60;

                success_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==9),3)/60; % stops
                failure_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==10),3)/60;
            elseif abc < 9
                success_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==7)/3,2)/60; % plate crosses
                failure_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==4)/3,2)/60;

                success_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==7)/3,3)/60; % stops
                failure_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==4)/3,3)/60;
            end
            
 
            success_markers_times(:,1) = success_markers_times(:,2)-2;  % first two seconds
            failure_markers_times(:,1) = failure_markers_times(:,2)-2;
            
            success_markers_times(:,6) = success_markers_times(:,4)+2;  % last two seconds
            failure_markers_times(:,6) = failure_markers_times(:,4)+2;
            
           success_idx = find(spiral_events(:,3)==7)/3;
            failure_idx = find(spiral_events(:,3)==4)/3;
            
            
              d = fields(avwt);
             avwt_all = getfield(avwt,d{1});
%              correct = find(diff(avwt_all')==0);
%              incorrect = find(50> abs(diff(avwt_all')) & abs(diff(avwt_all'))>0);
%              
%              succBad = find(isnan(sum(success_markers_times')));
%              failBad = find(isnan(sum(failure_markers_times')));
%              
%              idx = [find(isnan(tsCrossStop(:,2)))',correct(succBad), incorrect(failBad)];
%              avwt_all(unique(idx),:)=[];
            success_idx(isnan(sum(success_markers_times'))) = [];
            failure_idx(isnan(sum(failure_markers_times'))) = [];
            
            avwt_succ(isnan(sum(success_markers_times')),:) = [];
            avwt_fail(isnan(sum(failure_markers_times')),:) = [];
         
            success_markers_times(isnan(sum(success_markers_times')),:) = [];
            failure_markers_times(isnan(sum(failure_markers_times')),:) = [];
        
            f = find(diff(avwt_succ')~=0);
            ff = find(diff(avwt_fail')==0);
            avwt_succ(f,:)=[];
            success_idx(f)=[];
            success_markers_times(f,:)=[];
            avwt_fail(ff,:)=[];
            failure_idx(ff)=[];
            failure_markers_times(ff,:)=[];
            
            arena_ratemapper_lightscoring
%             arena_ratemapper
            save(['E:\nitzlabdata\new_ratemaps_12262012/parietal/' rats(abc,:) '_rec' recordings(r).name '_ratemap_anal_parietal.mat'])
            clear tfile* dvt* pos *spiral* avwt* success* failure*
        
            end
              end
% %         elseif isempty(hd)
% %         
        end

        
%         yourStr = '_parietal';
%         allVars = who;
%         for i = 1:length(allVars)
%         eval( [ allVars{i}, yourStr, ' = ', allVars{i} ';'] );
%         end
%         clearvars -except *parietal
        
%         save(['C:/Users/Nitz_Lab/David/new_ratemaps/parietal/' rats(abc).name '_rec' recordings(r).name '_ratemap_anal_parietal.mat'])
        end
        clear tfile* dvt* pos *spiral* *markers*
        
 end
    list_rats{abc} = list_to_ratemap;
    list_rats_names{abc} = list_to_ratemap_name;
    list_need{abc} = list_need_files;
end

 
