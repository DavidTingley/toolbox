clear 
rats = ['AA8';'DN3';'NS1'; 'NS2';'NS3';'NS4';'NS5';'NS8'];
% unused = dir('F:\David\*');
% med_only = dir('F:\David\*');
p = 'E:\nitzlabdata';
unused_recording = 0;
% medium_recording = 0;
% ratwire_loc_list = [zeros(1,16),17:48; ...  % AA8
%     1:16 ,zeros(1,16),33:48; ...           % DN3
%     zeros(1,16),zeros(1,16),zeros(1,16); ... % NS1
%     zeros(1,16),17:32,zeros(1,16); ...      % NS2
%     zeros(1,16),17:32,zeros(1,16); ...      % NS3
%     zeros(1,16),17:32,zeros(1,16); ...      % NS4
%     zeros(1,16),17:48;                    % NS5
%     zeros(1,16),17:48];                    % NS8
                    % this lists the wires that ARE NOT in the BF for each rat
                    % order is AA8 DN3 NS1 NS2 NS3 NS4
    ratwire_loc_list = [zeros(1,16),17:32,zeros(1,16); ...  % AA8
    1:16 ,zeros(1,16),33:48; ...           % DN3
    zeros(1,16),zeros(1,16),zeros(1,16); ... % NS1
    zeros(1,16),17:32,zeros(1,16); ...      % NS2
    zeros(1,16),17:32,zeros(1,16); ...      % NS3
    zeros(1,16),17:32,zeros(1,16); ...      % NS4
    zeros(1,16),17:48;                    % NS5
    zeros(1,16),17:48];                    % NS8

%    for iteration = 1:20                  

% getting rats directory
for abc = 1:8
    ratpath = [p '/' sprintf(rats(abc,:)) '/'];
    
    % getting recordings directory
    recordings = dir([ratpath '/*']);
    recordings = recordings(3:length(recordings));
    
     
    % recording check and make list of folders that need files
    [list_to_ratemap,list_to_ratemap_name, list_need_files] = rate_map_check(ratpath,recordings);
%     % makes a txt file that lists what is needed for each recording 
%     if abc == 1
%          f = fopen('list_needed_files.txt','w');
%     elseif abc > 1
%          f = fopen('list_needed_files.txt','a');
%     end
%      for i = 1:length(list_need_files)
%          s = []; 
%          if ~isempty(list_need_files{i})
%          for j = 1:3
%              if ~isempty(list_need_files{i}{j})
%              s = [list_need_files{i}{j} ', ' s];
%              end
%          end 
%          if ~isempty(s)
%      fprintf(f,[rats(abc,:) '_' recordings(i).name ' needs ' s '\n']);
%          end
%          end
%      end
%      fclose(f) 
 
            
    % ratemap recordings with all files
    for rr = 1:length(list_to_ratemap) 
        r = list_to_ratemap{rr};
        
        % do not ratemap \unused folder of previously ratemapped files
%         name = ([rats(abc,:) '_rec' recordings(r).name '_ratemap_anal.mat']);
        
%         for i = 3:length(unused)
%             if sum(unused(i).name == name) == length(name)
%                 if unused(i).name == name
%                 unused_recording = 1;
%                 end
%             end
%         end
%        pos = d(:,2:6);
%          f = fopen([ratpath recordings(r).name '\tfile_list.txt'],'w');
         
        [tfile_list]=textread([ratpath recordings(r).name '/tfile_list.txt'],'%s');
        
        %%%% have t   for i = 3:length(med_only)
%             if sum(med_only(i).name == name) == length(name)
%                 medium_recording = 1;
%             end
%         end
        if unused_recording == 0
        % setup variables for ratemapping
        spiral_info = dir([ratpath recordings(r).name '/' '*spiral*']);
        dvt = dir([ratpath recordings(r).name '/*1.DVT']);
        dvt2 = dir([ratpath recordings(r).name '/*2.DVT']);

   
         
        connect_issue = dir([ratpath recordings(r).name '/*reverse*']);
        if ~isempty(connect_issue)
            disp('connectors may be reversed')
            disp([ratpath recordings(r).name])
            if 'E:\nitzlabdata/NS3/38' ~= [ratpath recordings(r).name]
            return
            end
            ratwire_loc_list(5,:) = [zeros(1,length(1:32)),33:48];
        end
        if isempty(connect_issue)
            % reset ratwire_loc_list variable
            ratwire_loc_list(5,:) = [zeros(1,16),17:32,zeros(1,16)];
        end
         
        % Head direction
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
          

        load([ratpath recordings(r).name '/' rats(abc,:) '_rec' recordings(r).name...
            '_spiral_info.mat'])
        load([ratpath recordings(r).name '/' rats(abc,:) '_rec' recordings(r).name...
            '_platecross_and_stoptimes.mat'])

        pos = d(:,2:6);
%          f = fopen([ratpath recordings(r).name '\tfile_list.txt'],'w');
         
        [tfile_list]=textread([ratpath recordings(r).name '/tfile_list.txt'],'%s');
        
        %%%% have to pull out parietal wires here
        
        % the following code finds sig files that match parietal wires from
        % the rat_wire_loc_list variable using the tfile list and removes
        % them from the tfile_list variable, so that they are not
        % ratemapped.  Iventually they will be ratemapped seperately for
        % parietal analysis
        count = 1;
        for i = 1:length(tfile_list)
           for j = 1:length(ratwire_loc_list(abc,:))
               if ratwire_loc_list(abc,j) ~= 0
                 wire = findstr(tfile_list{i},num2str(ratwire_loc_list(abc,j),'%03d'));
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
        
        %%% END tfile_list editing for only BF cells
        if isempty(tfile_list) == 0  % if only cells were in different brain region,
                                     % skip ratemapping
            

        tfile_list_char(1:length(tfile_list),:)=char(tfile_list);

        
        spr = distInfo(rats(abc,:),recordings(r).name);
        datapath = [ratpath recordings(r).name '/'];
        
        t = dir([datapath '*.mat']);
        avwt = [];
        avwt_succ=[];
        for i = 1:length(t)
            if ~isempty(regexpi(t(i).name,'.*avwt.*'))
                avwt = load([datapath t(i).name]);
            end
            if ~isempty(regexpi(t(i).name,'.*error_anal.*'))
                avwt = load([datapath t(i).name],t(i).name(1:end-15));
            end
        end
      
        
        
        % ratemaps
        datapath = [ratpath recordings(r).name '/']
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
        if ~isempty(hd)
            % should be markers 
            list = who; 
            for i = 1:length(list)
            if length(list{i}) == length('success_beg')
            if list{i} == 'success_beg'
            s = eval(list{i});
            end
            end
            if length(list{i}) == length('succ_beg')
            if list{i} == 'succ_beg'
            s = eval(list{i});
            end
            end
            if length(list{i}) == length('succ_begs')
            if list{i} == 'succ_begs'
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
            if length(list{i}) == length('failure_beg')
            if list{i} == 'failure_beg'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('fail_beg')
            if list{i} == 'fail_beg'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('fail_begs')
            if list{i} == 'fail_begs'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('fail_begin')
            if list{i} == 'fail_begin'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('failure_begin')
            if list{i} == 'failure_begin'
            f = eval(list{i});
            end
            end
            if length(list{i}) == length('failure_beg')
            if list{i} == 'failure_beg'
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
            if length(list{i}) == length('succ_ends')
            if list{i} == 'succ_ends'
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
            if length(list{i}) == length('fail_ends')
            if list{i} == 'fail_ends'
            f_end = eval(list{i});
            end
            end
            if length(list{i}) == length('failure_end')
            if list{i} == 'failure_end'
            f_end = eval(list{i}); 
            end
            end
            if length(list{i}) == length('incorrect_end')
            if list{i} == 'incorrect_end'
            f_end = eval(list{i});
            end
            end
            end
            
            clear *arker* *end *beg* *DVT* *nogo*
            if r == 24 && abc == 4
                % remove 40 and 50 from spiral events
                spiral_events(199:201,3) = [1 8 6];
                spiral_events(244:246,3) = [1 8 6];
                f(24,:)=[];
            end
            if r == 8 && abc==5
                spiral_events(31:33,3) = [1 8 6];
                spiral_events(208:210,3) = [1 8 6];
                spiral_events(211:213,3) = [1 8 6];
                spiral_events(259:261,3) = [1 8 6];
            end
            if r == 11 && abc==5
                spiral_events(52:54,3) = [1 8 6]; 
                spiral_events(76:78,3) = [1 8 6];
                tsCrossStop = tsCrossStop+(46.7+1.413)*60;
            end
            if r == 9 && abc == 5
               tsCrossStop(11:end,2:3) = tsCrossStop(11:end,2:3)+13170;
            end
            if size(s,1) ~= length(spiral_events(find(spiral_events(:,3)==9),2))...
                    && size(f,1) ~= length(spiral_events(find(spiral_events(:,3)==10),2))
                error('something isnt matching')
            end
            
            if size(s,1) == length(spiral_events(find(spiral_events(:,3)==9),2))...
                    && size(f,1) == length(spiral_events(find(spiral_events(:,3)==10),2))
                disp(datapath)

            success_markers_times(:,2) = s(:,4);  % light flashes
            failure_markers_times(:,2) = f(:,4);
            
            success_markers_times(:,3) = spiral_events(find(spiral_events(:,3)==9),2); % nose-pokes
            failure_markers_times(:,3) = spiral_events(find(spiral_events(:,3)==10),2);
            
            success_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==7)/3,2)/60; % plate crosses
            failure_markers_times(:,4) = tsCrossStop(find(spiral_events(:,3)==4)/3,2)/60;
            
            success_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==7)/3,3)/60; % stops
            failure_markers_times(:,5) = tsCrossStop(find(spiral_events(:,3)==4)/3,3)/60;
 
            success_markers_times(:,1) = success_markers_times(:,2)-1;  % first two seconds
            failure_markers_times(:,1) = failure_markers_times(:,2)-1;
            
            success_markers_times(:,6) = success_markers_times(:,5)+1;  % last two seconds
            failure_markers_times(:,6) = failure_markers_times(:,5)+1;
            
            success_idx = find(spiral_events(:,3)==7)/3;
            failure_idx = find(spiral_events(:,3)==4)/3;
            
            
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
            end
   
            
            
        % checks for light-flash\trial missmatches 
        % if a mismatch is found this replaces the light flash with the
        % take off time
        if find(diff(success_markers_times')<0)
        for i = 1:size(success_markers_times,1)
           for j = 1:5
              if (success_markers_times(i,j+1)-success_markers_times(i,j)) <0
                 if j == 2
                    loc = find(spiral_events(:,2)==success_markers_times(i,3));
                    success_markers_times(i,j) = spiral_events(loc-1,2); 
                    success_markers_times(i,1) = spiral_events(loc-1,2)-2; 
                 end
              end
           end
        end
        end
        
        if find(diff(success_markers_times')>10)
        for i = 1:size(success_markers_times,1)
           for j = 1:5
              if (success_markers_times(i,j+1)-success_markers_times(i,j)) >10
                 if j == 2                    
                    loc = find(spiral_events(:,2)==success_markers_times(i,3));
                    success_markers_times(i,j) = spiral_events(loc-1,2); 
                    success_markers_times(i,1) = spiral_events(loc-1,2)-2; 
                    disp('change made')
                 end
              end
           end
        end
        end
                    %removes NaN's \ No Goes 
         if ~isempty(avwt)
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

        if size(avwt_succ,1) + size(avwt_fail,1) ~= size(success_markers_times,1) + ...
                size(failure_markers_times,1) & ~isempty(avwt)
            
            error('wtf')
        end
       
        if ~isempty(find(diff(success_markers_times')<0))
            error('time stamp mismatch')
        end

                
        TO_times{abc,rr}= spiral_events(find(spiral_events(:,3)==5),2);
        LF_times{abc,rr} =  s(:,4);
            % check for and remove NaNs that may have been introduced from
            % tsCrossStop times
%             temp = [];
%             for i = 1:size(success_markers_times,1)
%                 for j = 1:6
%                    if isnan(success_markers_times(i,j))
%                        temp = [temp; i];
%                    end
%                 end
%             end
%             success_markers_times(temp,:) = [];
%             temp = [];
%             for i = 1:size(failure_markers_times,1)
%                 for j = 1:6
%                    if isnan(failure_markers_times(i,j))
%                        temp = [temp;i];
%                    end
%                 end
%             end
%             failure_markers_times(temp,:) = [];
%             for temp = 1:size(success_markers_times,1)
%                 top = round(max(success_markers_times(:))-success_markers_times(temp,end));
%                 bot = round(min(success_markers_times(:))-success_markers_times(temp,1));
%                 success_markers_times(temp,:) = success_markers_times(temp,:) + randi([bot,top],1,1);
%             end
%             for temp = 1:size(failure_markers_times,1)
%                 top = round(max(failure_markers_times(:))-failure_markers_times(temp,end));
%                 bot = round(min(failure_markers_times(:))-failure_markers_times(temp,1));
%                 failure_markers_times(temp,:) = failure_markers_times(temp,:) + randi([bot,top],1,1);
%             end
            
            if strcmp([rats(abc,:) '_rec' recordings(r).name '_ratemap_anal'], 'NS2_rec21_ratemap_anal')
                disp('here we go');
            end
            arena_ratemapper_lightscoring
%         arena_ratemapper 

%         yourStr = '_parietal';
%         allVars = who;
%         for i = 1:length(allVars)
%         eval( [ allVars{i}, yourStr, ' = ', allVars{i} ';'] );
%         end
%         clearvars -except *parietal

% 
        save(['E:\nitzlabdata/new_ratemaps_12262012/' rats(abc,:) '_rec' recordings(r).name '_ratemap_anal.mat'])
%         if medium_recording == 1
%            save(['C:\Users\Nitz_Lab\David\new_ratemaps\medium_only\' rats(abc,:) '_rec' recordings(r).name '_ratemap_anal.mat'])
%         end      
        end
            end
        end
        clear tfile* dvt* pos *spiral* avwt* success* failure*
        
    end
    list_rats{abc} = list_to_ratemap;
    list_rats_names{abc} = list_to_ratemap_name;
    list_need{abc} = list_need_files;
    unused_recording = 0;
        end

    end
    end

%     end
 
