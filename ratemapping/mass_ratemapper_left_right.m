clear 
rats = dir('C:\Users\Nitz_Lab\David\*');
unused = dir('C:\Users\Nitz_Lab\David\new_ratemaps\unused\*');
p = pwd;
rats = rats(3:length(rats));
unused_recording = 0;
%% left brain
ratwire_loc_list_left = [zeros(1,16),17:48; ...  % AA8
    1:16 ,zeros(1,16),33:48; ...           % DN3
    1:16,zeros(1,16),33:48; ...          % NS1
    zeros(1,16),17:48; ...                  % NS2
    zeros(1,16),17:48; ...                  % NS3
    zeros(1,16),17:48; ...                  % NS4
    1:48;                                   % NS5
    1:48];                                  % NS8
                    % this lists the wires that ARE NOT in the LEFT BF for each rat
                    % order is AA8 DN3 NS1 NS2 NS3 NS4

% getting rats directory
for abc = [1:8]
    ratpath = [p '\' sprintf(rats(abc).name) '\'];
    
    % getting recordings directory
    recordings = dir([ratpath '\*']);
    recordings = recordings(3:length(recordings));
    
    
    % recording check and make list of folders that need files
    [list_to_ratemap,list_to_ratemap_name, list_need_files] = rate_map_check(ratpath,recordings);

    
    % ratemap recordings with all files
    for rr = 1:length(list_to_ratemap)
        r = list_to_ratemap{rr};
        spiral_info = dir([ratpath recordings(r).name '\' '*spiral*']);
        dvt = dir([ratpath recordings(r).name '\*1.DVT']);
        dvt2 = dir([ratpath recordings(r).name '\*2.DVT']);

        spr = distInfo(rats(abc).name,recordings(r).name);
        %do not ratemap unused
        name = ([rats(abc).name '_rec' recordings(r).name '_ratemap_anal.mat']);
        for i = 3:length(unused)
            if sum(unused(i).name == name) == length(name)
                if unused(i).name == name
                unused_recording = 1;
                end
            end
        end
        if unused_recording == 0
        connect_issue = dir([ratpath recordings(r).name '\*reverse*']);
        if ~isempty(connect_issue)
            disp('connectors may be reversed')
            disp([ratpath recordings(r).name])
            if 'C:\Users\Nitz_Lab\David\NS3\38' ~= [ratpath recordings(r).name]
            return
            end
            ratwire_loc_list(5,:) = [zeros(1,length(1:32)),33:48];
        end
        if isempty(connect_issue)
            % reset ratwire_loc_list variable
            ratwire_loc_list(5,:) = [zeros(1,16),17:32,zeros(1,16)];
        end
         
        % Head direction
                d = load([ratpath recordings(r).name '\' dvt.name]);
        if size(d,2) == 6
              %[headDir] = dvtHeadDir([ratpath recordings(r).name '\' dvt.name]);
              [headDir] = getHeadDirs([ratpath recordings(r).name '\' dvt.name]); 
        end
        if size(d,2) < 6
            d = load([ratpath recordings(r).name '\' dvt2.name]);
            %[headDir] = dvtHeadDir([ratpath recordings(r).name '\' dvt2.name]);
            [headDir] = getHeadDirs([ratpath recordings(r).name '\' dvt.name]); 
        end    
          

        load([ratpath recordings(r).name '\' rats(abc).name '_rec' recordings(r).name...
            '_spiral_info.mat'])

        pos = d(:,2:6);
        [tfile_list]=textread([ratpath recordings(r).name '\tfile_list.txt'],'%s');
        
        %%%% have to pull out parietal wires here
        
        % the following code finds sig files that match parietal wires from
        % the rat_wire_loc_list variable using the tfile list and removes
        % them from the tfile_list variable, so that they are not
        % ratemapped.  Eventually they will be ratemapped seperately for
        % parietal analysis
        count = 1;
        for i = 1:length(tfile_list)
           for j = 1:length(ratwire_loc_list_left(abc,:))
               if ratwire_loc_list_left(abc,j) ~= 0
                 wire = findstr(tfile_list{i},num2str(ratwire_loc_list_left(abc,j),'%03d'));
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
        if isempty(tfile_list) == 0  % if only cells were parietal, skip ratemapping
        tfile_list_char(1:length(tfile_list),:)=char(tfile_list);

        datapath = [ratpath recordings(r).name '\'];
        
        t = dir([datapath '*.mat']);
        for i = 1:length(t)
            if ~isempty(regexpi(t(i).name,'.*avwt.*'))
                avwt = load([datapath t(i).name]);
            end
        end
        
        
        arena_ratemapper 
        
%         yourStr = '_parietal';
%         allVars = who;
%         for i = 1:length(allVars)
%         eval( [ allVars{i}, yourStr, ' = ', allVars{i} ';'] );
%         end
%         clearvars -except *parietal
        
        save(['C:\Users\Nitz_Lab\David\new_ratemaps\left\' rats(abc).name '_rec' recordings(r).name '_ratemap_anal_left.mat'])
        end
        clear tfile* dvt* pos *spiral*
        end

    list_rats{abc} = list_to_ratemap;
    list_rats_names{abc} = list_to_ratemap_name;
    list_need{abc} = list_need_files;
    unused_recording = 0;
    end
end


%% right brain
 unused_recording = 0;
ratwire_loc_list_right = [1:48; ...         % AA8
    1:48; ...                               % DN3
    zeros(1,16),17:48; ...                  % NS1
    1:32,zeros(1,16); ...                   % NS2
    1:32,zeros(1,16); ...                   % NS3
    1:32,zeros(1,16); ...                   % NS4
    zeros(1,16),17:48;                      % NS5
    zeros(1,16),17:48;];                    % NS8
                    % this lists the wires that ARE NOT in the RIGHT BF for each rat
                    % order is AA8 DN3 NS1 NS2 NS3 NS4

% getting rats directory
for abc = [1:8]
    ratpath = [p '\' sprintf(rats(abc).name) '\'];
    
    % getting recordings directory
    recordings = dir([ratpath '\*']);
    recordings = recordings(3:length(recordings));
    
    
    % recording check and make list of folders that need files
    [list_to_ratemap,list_to_ratemap_name, list_need_files] = rate_map_check(ratpath,recordings);

    
    % ratemap recordings with all files
    for rr = 1:length(list_to_ratemap)
        r = list_to_ratemap{rr};
        spiral_info = dir([ratpath recordings(r).name '\' '*spiral*']);
        dvt = dir([ratpath recordings(r).name '\*1.DVT']);
        dvt2 = dir([ratpath recordings(r).name '\*2.DVT']);

        spr = distInfo(rats(abc).name,recordings(r).name);
        % do not ratemap unused
        name = ([rats(abc).name '_rec' recordings(r).name '_ratemap_anal.mat']);
        for i = 3:length(unused)
            if sum(unused(i).name == name) == length(name)
                if unused(i).name == name
                unused_recording = 1;
                end
            end
        end
        if unused_recording == 0
            
        connect_issue = dir([ratpath recordings(r).name '\*reverse*']);
        if ~isempty(connect_issue)
            disp('connectors may be reversed')
            disp([ratpath recordings(r).name])
            if 'C:\Users\Nitz_Lab\David\NS3\38' ~= [ratpath recordings(r).name]
            return
            end
            ratwire_loc_list(5,:) = [zeros(1,length(1:32)),33:48];
        end
        if isempty(connect_issue)
            % reset ratwire_loc_list variable
            ratwire_loc_list(5,:) = [zeros(1,16),17:32,zeros(1,16)];
        end
         
        % Head direction
                d = load([ratpath recordings(r).name '\' dvt.name]);
        if size(d,2) == 6
              [headDir] = dvtHeadDir([ratpath recordings(r).name '\' dvt.name]);
        end
        if size(d,2) < 6
            d = load([ratpath recordings(r).name '\' dvt2.name]);
            [headDir] = dvtHeadDir([ratpath recordings(r).name '\' dvt2.name]);
        end    
          

        load([ratpath recordings(r).name '\' rats(abc).name '_rec' recordings(r).name...
            '_spiral_info.mat'])

        pos = d(:,2:6);
        [tfile_list]=textread([ratpath recordings(r).name '\tfile_list.txt'],'%s');
        
        %%%% have to pull out parietal wires here
        
        % the following code finds sig files that match parietal wires from
        % the rat_wire_loc_list variable using the tfile list and removes
        % them from the tfile_list variable, so that they are not
        % ratemapped.  Eventually they will be ratemapped seperately for
        % parietal analysis
        count = 1;
        for i = 1:length(tfile_list)
           for j = 1:length(ratwire_loc_list_right(abc,:))
               if ratwire_loc_list_right(abc,j) ~= 0
                 wire = findstr(tfile_list{i},num2str(ratwire_loc_list_right(abc,j),'%03d'));
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

        tfile_list_char(1:length(tfile_list),:)=char(tfile_list);
        if isempty(tfile_list) == 0  % if only cells were parietal, skip ratemapping
        datapath = [ratpath recordings(r).name '\'];
        
        t = dir([datapath '*.mat']);
        for i = 1:length(t)
            if ~isempty(regexpi(t(i).name,'.*avwt.*'))
                avwt = load([datapath t(i).name]);
            end
        end
        
        
        arena_ratemapper 
        
%         yourStr = '_parietal';
%         allVars = who;
%         for i = 1:length(allVars)
%         eval( [ allVars{i}, yourStr, ' = ', allVars{i} ';'] );
%         end
%         clearvars -except *parietal
        
        save(['C:\Users\Nitz_Lab\David\new_ratemaps\right\' rats(abc).name '_rec' recordings(r).name '_ratemap_anal_right.mat'])
        end
        clear tfile* dvt* pos *spiral*
        end

    list_rats{abc} = list_to_ratemap;
    list_rats_names{abc} = list_to_ratemap_name;
    list_need{abc} = list_need_files;
    unused_recording = 0;
    end
end