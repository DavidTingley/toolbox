

dvt =  AG1_rec01_030113001_2;
right = load('Event003');
left = load('Event004');
tfile_list = dir('sig0*');



% scatter(dvt(:,3),dvt(:,4)); axis([0 900 0 900])
% [x y] = ginput(4);

start_left_success_trial = [];
end_left_success_trial = [];
start_left_failure_trial = [];
end_left_failure_trial = [];

start_right_success_trial = [];
end_right_success_trial = [];
start_right_failure_trial = [];
end_right_failure_trial = [];

for i = 1:length(left)
    
    [ind time] = searchclosest(dvt(:,2),left(i))
    scatter(dvt(ind-300:ind+500,3),dvt(ind-300:ind+500,4)); axis([0 900 0 900])
    title('left trial')
      [leftx lefty button] = ginput(2);

      
      if button == 1
      [temp start_trial]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(1),lefty(1)]));
      start_left_success_trial = [start_left_success_trial;start_trial+ind-300];
      [temp end_trial ]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(2),lefty(2)]));
      end_left_success_trial = [end_left_success_trial;end_trial+ind-300];
      end
      
      if button == 3
      [temp start_trial]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(1),lefty(1)]));
      start_left_failure_trial = [start_left_failure_trial;start_trial+ind-300];
      [temp end_trial ]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(2),lefty(2)]));
      end_left_failure_trial = [end_left_failure_trial;end_trial+ind-300];
      end
end


for i = 1:length(right)  
    [ind time] = searchclosest(dvt(:,2),right(i))
    scatter(dvt(ind-300:ind+500,3),dvt(ind-300:ind+500,4)); axis([0 900 0 900])
    title('right trial')
      [rightx righty button] = ginput(2);  
    
      if button == 1
      [temp start_trial]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(1),lefty(1)]));
      start_right_success_trial = [start_right_success_trial;start_trial+ind-300];
      [temp end_trial ]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(2),lefty(2)]));
      end_right_success_trial = [end_right_success_trial;end_trial+ind-300];
      end
      
      if button == 3
      [temp start_trial]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(1),lefty(1)]));
      start_right_failure_trial = [start_right_failure_trial;start_trial+ind-300];
      [temp end_trial ]= min(distance(dvt(ind-300:ind+500,3:4),[leftx(2),lefty(2)]));
      end_right_failure_trial = [end_right_failure_trial;end_trial+ind-300];
      end
end



