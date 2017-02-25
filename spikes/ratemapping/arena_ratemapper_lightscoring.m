
%2D ratemap smoothing function **************************************
half_narrow_val=7;
narrow_stdev_val=3.5;
half_narrow=half_narrow_val;
narrow_stdev=narrow_stdev_val;
[x y]=meshgrid(-half_narrow:1:half_narrow);
narrow_gaussian=double_gaussian(x,y,narrow_stdev);
narrow_gaussian=narrow_gaussian./sum(sum(narrow_gaussian));

%linear ratemap smoothing function *********************************
z=0:1:26;
y3=gaussmf(z,[2 12]);
y3=y3/sum(y3);

%peth smoothing function

z_peth=0:1:6;
y3_peth=gaussmf(z_peth,[1 3]);
y3_peth=y3_peth./sum(y3_peth);

%parsing behavioral scoring

no_cells=length(tfile_list(:,1));
% ratemaps=zeros(no_cells,1014,1014,2);

% add light flash time stamps
% success_starts=find(spiral_events(:,3)==5);
% failure_ends=find(spiral_events(:,3)==2);
% success_nose=find(spiral_events(:,3)==9);
% success_ends=find(spiral_events(:,3)==7);
% failure_nose=find(spiral_events(:,3)==10);
% failure_ends=find(spiral_events(:,3)==4);
% 
no_successes=length(success_markers_times);
no_failures=length(failure_markers_times);
% 
g=length(success_markers_times);
gg=length(success_markers_times);
ggg=length(success_markers_times);
h=length(failure_markers_times);
hh=length(failure_markers_times);
hhh=length(failure_markers_times);
if g==gg & gg==ggg
    success_markers_match=1
else
    potential_success_marker_match_problem=1
end
if h==hh & hh==hhh
    failure_markers_match=1
else
    potential_failure_marker_match_problem=1
end
% 
% success_markers_points=zeros(g,5);
% failure_markers_points=zeros(h,5);
% success_markers_times=zeros(g,5);
% failure_markers_times=zeros(h,5);
% for i=1:g
%     success_markers_points(i,:)=[spiral_events(success_markers_times(i),1)-120 spiral_events(success_markers_times(i),1) spiral_events(success_nose(i),1) spiral_events(success_ends(i),1) spiral_events(success_ends(i),1)+120];
% end
% for i=1:h
%     failure_markers_points(i,:)=[spiral_events(failure_markers_times(i),1)-120 spiral_events(failure_markers_times(i),1) spiral_events(failure_nose(i),1) spiral_events(failure_ends(i),1) spiral_events(failure_ends(i),1)+120];
% end
% for i=1:g
%     success_markers_times(i,:)=[spiral_events(success_markers_times(i),2)-2 spiral_events(success_markers_times(i),2) spiral_events(success_nose(i),2) spiral_events(success_ends(i),2) spiral_events(success_ends(i),2)+2];
% end
% for i=1:h
%     failure_markers_times(i,:)=[spiral_events(failure_markers_times(i),2)-2 spiral_events(failure_markers_times(i),2) spiral_events(failure_nose(i),2) spiral_events(failure_ends(i),2) spiral_events(failure_ends(i),2)+2];
% end
success_start_nose_times=success_markers_times(:,3)-success_markers_times(:,2);
success_nose_end_times=success_markers_times(:,4)-success_markers_times(:,3);
failure_start_nose_times=failure_markers_times(:,3)-failure_markers_times(:,2);
failure_nose_end_times=failure_markers_times(:,4)-failure_markers_times(:,3);

b=find(success_start_nose_times<0 | success_start_nose_times>5);
bb=find(failure_start_nose_times<0 | failure_start_nose_times>5);
bbb=find(success_nose_end_times<0 | success_nose_end_times>5);
bbbb=find(failure_nose_end_times<0 | failure_nose_end_times>5);
if length(b) >0
    problem_with_event_ordering=1
end
if length(bb) >0
    problem_with_event_ordering=1
end
if length(bbb) >0
    problem_with_event_ordering=1
end
if length(bbbb) >0
    problem_with_event_ordering=1
end

%filling in missing pos points
pos_ascii=pos;

for i=6:length(pos_ascii(:,1))-5
if pos_ascii(i,2)==1 & pos_ascii(i,3)==1
x=pos_ascii(i-5:i+5,2:3);
y=find(x(:,1)>1 & x(:,2)>1);
if length(y)>0
z=(y(:)-6);
zz=min(abs(z));
zzz=find(abs(y(:)-6)==zz);
if z(zzz(1))<0
xx=pos_ascii(i+z(zzz(1))+1:i+z(zzz(1))+5,2:3);
yy=find(xx(:,1)>1 & xx(:,2)>1);
x1=pos_ascii(i+z(zzz(1)),2);
y1=pos_ascii(i+z(zzz(1)),3);
if length(yy)>0
x2=pos_ascii(i+z(zzz(1))+1+yy(1),2);
y2=pos_ascii(i+z(zzz(1))+1+yy(1),3);
else
x2=pos_ascii(i+z(zzz(1)),2);
y2=pos_ascii(i+z(zzz(1)),3);
end
if abs(x1-x2)<15 & abs(y1-y2)<15
x3=round((x1+x2)/2);
y3=round((y1+y2)/2);
else
x3=x1;
y3=y1;
end
pos_ascii(i,2)=x3;
pos_ascii(i,3)=y3;
elseif z(zzz(1))>0
xx=pos_ascii(i+z(zzz(1))-6:i+z(zzz(1))-1,2:3);
yy=find(xx(:,1)>1 & xx(:,2)>1);
x1=pos_ascii(i+z(zzz(1)),2);
y1=pos_ascii(i+z(zzz(1)),3);
if length(yy)>0
x2=pos_ascii(i+z(zzz(1))-6+xx(length(xx)),2);
y2=pos_ascii(i+z(zzz(1))-6+xx(length(xx)),3);
else
x2=pos_ascii(i+z(zzz(1)),2);
y2=pos_ascii(i+z(zzz(1)),3);
end
if abs(x1-x2)<15 & abs(y1-y2)<15
x3=round((x1+x2)/2);
y3=round((y1+y2)/2);
else
x3=x1;
y3=y1;
end
pos_ascii(i,2)=x3;
pos_ascii(i,3)=y3;
end
end
end
end

done='filled_in_missing_pos_points'

pos=pos_ascii;

% calculating velocity and acceleration for each position point

vels=zeros(length(pos_ascii(:,1)),1);
for i=3:length(pos_ascii)-3
if pos_ascii(i+2,2)>1 & pos_ascii(i-2,2)>1 & pos_ascii(i+2,3)>1 & pos_ascii(i-2,3)>1
x_diff=abs(pos_ascii(i+2,2)-pos_ascii(i-2,2));
y_diff=abs(pos_ascii(i+2,3)-pos_ascii(i-2,3));
else
x_diff=1000;
y_diff=1000;
end
if x_diff==0 & y_diff==0
vels(i)=0;
elseif y_diff==0 & x_diff>0
vels(i)=x_diff;
elseif x_diff==0 & y_diff>0
vels(i)=y_diff;
else 
vels(i)=y_diff/sin(atan(y_diff/x_diff));
end
end
vels(1:4)=mean(vels(5:10));
vels(length(pos_ascii(:,1))-4:length(pos_ascii(:,1)))=mean(vels(length(pos_ascii(:,1))-10:length(pos_ascii(:,1))-5));
vels=round(vels);
vels=vels+1;
z=0:1:13;
y3=gaussmf(z,[1 6]);
y3=y3/sum(y3);
vels_filt=filtfilt(y3,1,vels);
x=diff(vels_filt);
accs=zeros(length(vels_filt),1);
accs(2:length(accs),1)=x;
accs(1,1)=accs(2,1);

done='vel_acc_calcs'

pos=pos_ascii;

% making 2D ratemaps

% for loop=1:no_cells
% 
% pos_spikes_success=zeros(1000,1000);
% pos_spikes_failure=zeros(1000,1000);
% pos_occ_success=zeros(1000,1000);
% pos_occ_failure=zeros(1000,1000);
% pos_rates_success=zeros(1000,1000);
% pos_rates_failure=zeros(1000,1000);
% 
% for i=1:length(success_markers_times)
% for j=spiral_events(success_markers_times(i),1)-120:spiral_events(success_ends(i),1)+120
% pos_occ_success(pos(j,2),pos(j,3))=pos_occ_success(pos(j,2),pos(j,3))+1;
% end
% end
% 
% for i=1:length(failure_markers_times)
% for j=spiral_events(failure_markers_times(i),1)-120:spiral_events(failure_ends(i),1)+120
% pos_occ_failure(pos(j,2),pos(j,3))=pos_occ_failure(pos(j,2),pos(j,3))+1;
% end
% end
% 
% tfile=load(tfile_list_char(loop,:));
% 
% for i=1:length(success_markers_times)
% for j=spiral_events(success_markers_times(i),1)-120:spiral_events(success_ends(i),1)+120
% m=find(tfile>=pos(j,1)-(.0167/2) & tfile<pos(j,1)+(.0167/2));
% mm=length(m);
% pos_spikes_success(pos(j,2),pos(j,3))=pos_spikes_success(pos(j,2),pos(j,3))+mm;
% end
% end
% 
% for i=1:length(failure_markers_times)
% for j=spiral_events(failure_markers_times(i),1)-120:spiral_events(failure_ends(i),1)+120
% m=find(tfile>=pos(j,1)-(.0167/2) & tfile<pos(j,1)+(.0167/2));
% mm=length(m);
% pos_spikes_failure(pos(j,2),pos(j,3))=pos_spikes_failure(pos(j,2),pos(j,3))+mm;
% end
% end
% 
% for i=1:1000
% for j=1:1000
% if pos_occ_success(i,j)>0
% pos_rates_success(i,j)=(pos_spikes_success(i,j)/pos_occ_success(i,j))*60;
% else
% pos_rates_success(i,j)=-.5;
% end
% if pos_occ_failure(i,j)>0
% pos_rates_failure(i,j)=(pos_spikes_failure(i,j)/pos_occ_failure(i,j))*60;
% else
% pos_rates_failure(i,j)=-.5;
% end
% end
% end
% 
% pos_rates_success=conv2(pos_rates_success,narrow_gaussian);
% pos_rates_failure=conv2(pos_rates_failure,narrow_gaussian);
% 
% ratemaps(loop,:,:,1)=pos_rates_success;
% ratemaps(loop,:,:,2)=pos_rates_failure;

% twoD_ratemap_completed=loop
% end

%making peri-event time histograms - 50ms bins going back and fwd 1 sec
%from start, nose, and end markers - note that bin 21 is 50ms period
%surrounding the marker itself

% success_markers_times_peth=zeros(no_cells,length(success_markers_times),41);
% failure_markers_times_peth=zeros(no_cells,length(failure_markers_times),41);
% success_nose_peth=zeros(no_cells,length(success_nose),41);
% failure_nose_peth=zeros(no_cells,length(failure_nose),41);
% success_ends_peth=zeros(no_cells,length(success_ends),41);
% failure_ends_peth=zeros(no_cells,length(failure_ends),41);
% success_markers_times_peth_smooth=zeros(no_cells,length(success_markers_times),41);
% failure_markers_times_peth_smooth=zeros(no_cells,length(failure_markers_times),41);
% success_nose_peth_smooth=zeros(no_cells,length(success_nose),41);
% failure_nose_peth_smooth=zeros(no_cells,length(failure_nose),41);
% success_ends_peth_smooth=zeros(no_cells,length(success_ends),41);
% failure_ends_peth_smooth=zeros(no_cells,length(failure_ends),41);
% 
% for loop=1:no_cells
% tfile=load([datapath tfile_list_char(loop,:)]);
% for i=1:length(success_markers_times)
%     for j=-20:20
%     x=find(tfile>(success_markers_times(i,2)+((j*.050)-.025)) & tfile<(success_markers_times(i,2)+((j*.050)+.0249999)));
%     success_markers_times_peth(loop,i,j+21)=length(x)*20;
%     x=find(tfile>(success_markers_times(i,3)+((j*.050)-.025)) & tfile<(success_markers_times(i,3)+((j*.050)+.0249999)));
%     success_nose_peth(loop,i,j+21)=length(x)*20;
%     x=find(tfile>(success_markers_times(i,4)+((j*.050)-.025)) & tfile<(success_markers_times(i,4)+((j*.050)+.0249999)));
%     success_ends_peth(loop,i,j+21)=length(x)*20;
%     end
% end
% for i=1:length(failure_markers_times)
%     for j=-20:20
%     x=find(tfile>(failure_markers_times(i,2)+((j*.050)-.025)) & tfile<(failure_markers_times(i,2)+((j*.050)+.0249999)));
%     failure_markers_times_peth(loop,i,j+21)=length(x)*20;
%     x=find(tfile>(failure_markers_times(i,3)+((j*.050)-.025)) & tfile<(failure_markers_times(i,3)+((j*.050)+.0249999)));
%     failure_nose_peth(loop,i,j+21)=length(x)*20;
%     x=find(tfile>(failure_markers_times(i,4)+((j*.050)-.025)) & tfile<(failure_markers_times(i,4)+((j*.050)+.0249999)));
%     failure_ends_peth(loop,i,j+21)=length(x)*20;
%     end
% end
% end
% 
% % making smoothed peths
% 
% for loop=1:no_cells
%     for j=1:no_successes
% success_markers_times_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,success_markers_times_peth(loop,j,:));
% success_nose_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,success_nose_peth(loop,j,:));
% success_ends_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,success_ends_peth(loop,j,:));
%     end
% end
% for loop=1:no_cells
%     for j=1:no_failures
% failure_markers_times_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,failure_markers_times_peth(loop,j,:));
% failure_nose_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,failure_nose_peth(loop,j,:));
% failure_ends_peth_smooth(loop,j,:)=filtfilt(y3_peth,1,failure_ends_peth(loop,j,:));
%     end
% end
% 
% % calc of means, stds, stes of peths
% 
% success_markers_times_peth_means_stds_stes=zeros(no_cells,41,3);
% failure_markers_times_peth_means_stds_stes=zeros(no_cells,41,3);
% success_nose_peth_means_stds_stes=zeros(no_cells,41,3);
% failure_nose_peth_means_stds_stes=zeros(no_cells,41,3);
% success_ends_peth_means_stds_stes=zeros(no_cells,41,3);
% failure_ends_peth_means_stds_stes=zeros(no_cells,41,3);
% 
% success_markers_times_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% failure_markers_times_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% success_nose_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% failure_nose_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% success_ends_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% failure_ends_peth_smooth_means_stds_stes=zeros(no_cells,41,3);
% 
% for loop=1:no_cells
%     for j=1:41
%     success_markers_times_peth_means_stds_stes(loop,j,1)=mean(success_markers_times_peth(loop,:,j));
%     success_markers_times_peth_means_stds_stes(loop,j,2)=std(success_markers_times_peth(loop,:,j));
%     success_markers_times_peth_means_stds_stes(loop,j,3)=std(success_markers_times_peth(loop,:,j))/(sqrt(length(success_markers_times)-1));
%     success_nose_peth_means_stds_stes(loop,j,1)=mean(success_nose_peth(loop,:,j));
%     success_nose_peth_means_stds_stes(loop,j,2)=std(success_nose_peth(loop,:,j));
%     success_nose_peth_means_stds_stes(loop,j,3)=std(success_nose_peth(loop,:,j))/(sqrt(length(success_nose)-1));
%     success_ends_peth_means_stds_stes(loop,j,1)=mean(success_ends_peth(loop,:,j));
%     success_ends_peth_means_stds_stes(loop,j,2)=std(success_ends_peth(loop,:,j));
%     success_ends_peth_means_stds_stes(loop,j,3)=std(success_ends_peth(loop,:,j))/(sqrt(length(success_ends)-1));
%     failure_markers_times_peth_means_stds_stes(loop,j,1)=mean(failure_markers_times_peth(loop,:,j));
%     failure_markers_times_peth_means_stds_stes(loop,j,2)=std(failure_markers_times_peth(loop,:,j));
%     failure_markers_times_peth_means_stds_stes(loop,j,3)=std(failure_markers_times_peth(loop,:,j))/(sqrt(length(failure_markers_times)-1));
%     failure_nose_peth_means_stds_stes(loop,j,1)=mean(failure_nose_peth(loop,:,j));
%     failure_nose_peth_means_stds_stes(loop,j,2)=std(failure_nose_peth(loop,:,j));
%     failure_nose_peth_means_stds_stes(loop,j,3)=std(failure_nose_peth(loop,:,j))/(sqrt(length(failure_nose)-1));
%     failure_ends_peth_means_stds_stes(loop,j,1)=mean(failure_ends_peth(loop,:,j));
%     failure_ends_peth_means_stds_stes(loop,j,2)=std(failure_ends_peth(loop,:,j));
%     failure_ends_peth_means_stds_stes(loop,j,3)=std(failure_ends_peth(loop,:,j))/(sqrt(length(failure_ends)-1));
%     
%     success_markers_times_peth_smooth_means_stds_stes(loop,j,1)=mean(success_markers_times_peth_smooth(loop,:,j));
%     success_markers_times_peth_smooth_means_stds_stes(loop,j,2)=std(success_markers_times_peth_smooth(loop,:,j));
%     success_markers_times_peth_smooth_means_stds_stes(loop,j,3)=std(success_markers_times_peth_smooth(loop,:,j))/(sqrt(length(success_markers_times_peth_smooth)-1));
%     success_nose_peth_smooth_means_stds_stes(loop,j,1)=mean(success_nose_peth_smooth(loop,:,j));
%     success_nose_peth_smooth_means_stds_stes(loop,j,3)=std(success_nose_peth_smooth(loop,:,j));
%     success_nose_peth_smooth_means_stds_stes(loop,j,3)=std(success_nose_peth_smooth(loop,:,j))/(sqrt(length(success_nose_peth_smooth)-1));
%     success_ends_peth_smooth_means_stds_stes(loop,j,1)=mean(success_ends_peth_smooth(loop,:,j));
%     success_ends_peth_smooth_means_stds_stes(loop,j,2)=std(success_ends_peth_smooth(loop,:,j));
%     success_ends_peth_smooth_means_stds_stes(loop,j,3)=std(success_ends_peth_smooth(loop,:,j))/(sqrt(length(success_ends_peth_smooth)-1));
%     failure_markers_times_peth_smooth_means_stds_stes(loop,j,1)=mean(failure_markers_times_peth_smooth(loop,:,j));
%     failure_markers_times_peth_smooth_means_stds_stes(loop,j,2)=std(failure_markers_times_peth_smooth(loop,:,j));
%     failure_markers_times_peth_smooth_means_stds_stes(loop,j,3)=std(failure_markers_times_peth_smooth(loop,:,j))/(sqrt(length(failure_markers_times_peth_smooth)-1));
%     failure_nose_smooth_peth_means_stds_stes(loop,j,1)=mean(failure_nose_peth_smooth(loop,:,j));
%     failure_nose_smooth_peth_means_stds_stes(loop,j,2)=std(failure_nose_peth_smooth(loop,:,j));
%     failure_nose_smooth_peth_means_stds_stes(loop,j,3)=std(failure_nose_peth_smooth(loop,:,j))/(sqrt(length(failure_nose_peth_smooth)-1));
%     failure_ends_smooth_peth_means_stds_stes(loop,j,1)=mean(failure_ends_peth_smooth(loop,:,j));
%     failure_ends_smooth_peth_means_stds_stes(loop,j,2)=std(failure_ends_peth_smooth(loop,:,j));
%     failure_ends_smooth_peth_means_stds_stes(loop,j,3)=std(failure_ends_peth_smooth(loop,:,j))/(sqrt(length(failure_ends_peth_smooth)-1));    
%     end
% end

%calc mean rate overall of each cell

mean_rates_overall=zeros(no_cells,1);
for loop=1:no_cells
    t = [datapath tfile_list_char(loop,:)];
    tfile=load(deblank(t));
    mean_rates_overall(loop)=length(tfile)/(tfile(length(tfile))-tfile(1));
end

%making a magnitude peth for each cell (a-b/a+b) with mean overall rate as
%the 'b' value

% success_markers_times_peth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% success_nose_peth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% success_ends_peth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% success_markers_times_peth_smooth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% success_nose_peth_smooth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% success_ends_peth_smooth_normed_mag=zeros(no_cells,length(success_markers_times),41);
% failure_markers_times_peth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% failure_nose_peth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% failure_ends_peth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% failure_markers_times_peth_smooth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% failure_nose_peth_smooth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% failure_ends_peth_smooth_normed_mag=zeros(no_cells,length(failure_markers_times),41);
% for loop=1:no_cells
%     success_markers_times_peth_normed_mag(loop,:,:)=(success_markers_times_peth(loop,:,:)-mean_rates_overall(loop))./(success_markers_times_peth(loop,:,:)+mean_rates_overall(loop));
%     success_nose_peth_normed_mag(loop,:,:)=(success_nose_peth(loop,:,:)-mean_rates_overall(loop))./(success_nose_peth(loop,:,:)+mean_rates_overall(loop));
%     success_ends_peth_normed_mag(loop,:,:)=(success_ends_peth(loop,:,:)-mean_rates_overall(loop))./(success_ends_peth(loop,:,:)+mean_rates_overall(loop));
%     failure_markers_times_peth_normed_mag(loop,:,:)=(failure_markers_times_peth(loop,:,:)-mean_rates_overall(loop))./(failure_markers_times_peth(loop,:,:)+mean_rates_overall(loop));
%     failure_nose_peth_normed_mag(loop,:,:)=(failure_nose_peth(loop,:,:)-mean_rates_overall(loop))./(failure_nose_peth(loop,:,:)+mean_rates_overall(loop));
%     failure_ends_peth_normed_mag(loop,:,:)=(failure_ends_peth(loop,:,:)-mean_rates_overall(loop))./(failure_ends_peth(loop,:,:)+mean_rates_overall(loop));
%     success_markers_times_peth_smooth_normed_mag(loop,:,:)=(success_markers_times_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(success_markers_times_peth_smooth(loop,:,:)+mean_rates_overall(loop));
%     success_nose_peth_smooth_normed_mag(loop,:,:)=(success_nose_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(success_nose_peth_smooth(loop,:,:)+mean_rates_overall(loop));
%     success_ends_peth_smooth_normed_mag(loop,:,:)=(success_ends_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(success_ends_peth_smooth(loop,:,:)+mean_rates_overall(loop));
%     failure_markers_times_peth_smooth_normed_mag(loop,:,:)=(failure_markers_times_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(failure_markers_times_peth_smooth(loop,:,:)+mean_rates_overall(loop));
%     failure_nose_peth_smooth_normed_mag(loop,:,:)=(failure_nose_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(failure_nose_peth_smooth(loop,:,:)+mean_rates_overall(loop));
%     failure_ends_peth_smooth_normed_mag(loop,:,:)=(failure_ends_peth_smooth(loop,:,:)-mean_rates_overall(loop))./(failure_ends_peth_smooth(loop,:,:)+mean_rates_overall(loop));
% end
    
%making a time-normalized 'ratemap' for the time preceding starts by 2 sec
%to the time following ends by 2 sec

success_timenorm_rates=zeros(no_cells,length(success_markers_times),54);
failure_timenorm_rates=zeros(no_cells,length(failure_markers_times),54);

for loop=1:no_cells
    tfile=load(deblank([datapath tfile_list_char(loop,:)]));
    for i=1:length(success_markers_times)
        % first 2 seconds
        for j=-12:-1
            x=find(tfile>(success_markers_times(i,2)+(j*.087)) & tfile<(success_markers_times(i,2)+(j*.087)+.087));
            success_timenorm_rates(loop,i,j+13)=length(x)*11.3636;
        end
        % stop to end
        for j=1:12
            x=find(tfile>(success_markers_times(i,5)+((j-1)*.087)) & tfile<(success_markers_times(i,5)+(((j-1)*.087)+.087)));
            success_timenorm_rates(loop,i,j+42)=length(x)*11.3636;
        end    
        w=success_markers_times(:,3)-success_markers_times(:,2);
        ww=success_markers_times(:,4)-success_markers_times(:,3);
        www=success_markers_times(:,5)-success_markers_times(:,4);
        run_out_times_succ = w;
        run_in_times_succ = ww;
        run_to_stop_times_succ = www;
        trial_timebin=w(i)/8;
        % run out
        for j=1:8
            x=find(tfile>success_markers_times(i,2)+((j-1)*trial_timebin) & tfile<(success_markers_times(i,2)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            success_timenorm_rates(loop,i,j+12)=length(x)*(1/trial_timebin);
        end
        % plate-cross
        trial_timebin=ww(i)/16;
        for j=1:16
            x=find(tfile>success_markers_times(i,3)+((j-1)*trial_timebin) & tfile<(success_markers_times(i,3)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            success_timenorm_rates(loop,i,j+20)=length(x)*(1/trial_timebin);
        end
        % stop
        trial_timebin=www(i)/6;
        for j=1:6
            x=find(tfile>success_markers_times(i,4)+((j-1)*trial_timebin) & tfile<(success_markers_times(i,4)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            success_timenorm_rates(loop,i,j+36)=length(x)*(1/trial_timebin);
        end
    end
end

for loop=1:no_cells
    tfile=load(deblank([datapath tfile_list_char(loop,:)]));
    for i=1:size(failure_markers_times,1)
        % first 2 seconds
        for j=-12:-1
            x=find(tfile>(failure_markers_times(i,2)+(j*.087)) & tfile<(failure_markers_times(i,2)+(j*.087)+.087));
            failure_timenorm_rates(loop,i,j+13)=length(x)*11.3636;
        end
       % stop to end
        for j=1:12
            x=find(tfile>(failure_markers_times(i,5)+((j-1)*.087)) & tfile<(failure_markers_times(i,5)+(((j-1)*.087)+.087)));
            failure_timenorm_rates(loop,i,j+42)=length(x)*11.3636;
        end    
        w=failure_markers_times(:,3)-failure_markers_times(:,2);
        ww=failure_markers_times(:,4)-failure_markers_times(:,3);
        www=failure_markers_times(:,5)-failure_markers_times(:,4);
        run_out_times_fail = w;
        run_in_times_fail = ww;
        run_to_stop_times_fail = www;
        trial_timebin=w(i)/8;
        % run out
        for j=1:8
            x=find(tfile>failure_markers_times(i,2)+((j-1)*trial_timebin) & tfile<(failure_markers_times(i,2)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            failure_timenorm_rates(loop,i,j+12)=length(x)*(1/trial_timebin);
        end
        % plate-cross
        trial_timebin=ww(i)/16;
        for j=1:16
            x=find(tfile>failure_markers_times(i,3)+((j-1)*trial_timebin) & tfile<(failure_markers_times(i,3)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            failure_timenorm_rates(loop,i,j+20)=length(x)*(1/trial_timebin);
        end
        % stop
        trial_timebin=www(i)/6;
        for j=1:6
            x=find(tfile>failure_markers_times(i,4)+((j-1)*trial_timebin) & tfile<(failure_markers_times(i,4)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
            failure_timenorm_rates(loop,i,j+36)=length(x)*(1/trial_timebin);
        end
    end
end

success_timenorm_rates_smooth=zeros(no_cells,length(success_markers_times),54);
failure_timenorm_rates_smooth=zeros(no_cells,length(failure_markers_times),54);

for loop=1:no_cells
    for i=1:no_successes
        success_timenorm_rates_smooth(loop,i,:)=filtfilt(y3_peth,1,success_timenorm_rates(loop,i,:));
    end
    for i=1:no_failures
        failure_timenorm_rates_smooth(loop,i,:)=filtfilt(y3_peth,1,failure_timenorm_rates(loop,i,:));
    end
end

%take means stds stes of smoothed timenorm rates

success_timenorm_rates_smooth_means_stds_stes=zeros(no_cells,54,3);
failure_timenorm_rates_smooth_means_stds_stes=zeros(no_cells,54,3);

for loop=1:no_cells
    for i=1:54
    success_timenorm_rates_smooth_means_stds_stes(loop,i,1)=mean(success_timenorm_rates_smooth(loop,:,i));
    success_timenorm_rates_smooth_means_stds_stes(loop,i,2)=std(success_timenorm_rates_smooth(loop,:,i));
    success_timenorm_rates_smooth_means_stds_stes(loop,i,3)=std(success_timenorm_rates_smooth(loop,:,i))/(sqrt(no_successes-1));
    failure_timenorm_rates_smooth_means_stds_stes(loop,i,1)=mean(failure_timenorm_rates_smooth(loop,:,i));
    failure_timenorm_rates_smooth_means_stds_stes(loop,i,2)=std(failure_timenorm_rates_smooth(loop,:,i));
    failure_timenorm_rates_smooth_means_stds_stes(loop,i,3)=std(failure_timenorm_rates_smooth(loop,:,i))/(sqrt(no_failures-1));
    end
end

%take diffs and mag of diffs between timenorm success and failure patterns

% diff_magdiffs_success_failure_mean_std_timenorm_rates=zeros(no_cells,54,4);
% 
% for loop=1:no_cells
%     diff_magdiffs_success_failure_mean_std_timenorm_rates(loop,:,1)=success_timenorm_rates_smooth_means_stds_stes(loop,:,1)-failure_timenorm_rates_smooth_means_stds_stes(loop,:,1);
%     diff_magdiffs_success_failure_mean_std_timenorm_rates(loop,:,2)=success_timenorm_rates_smooth_means_stds_stes(loop,:,2)-failure_timenorm_rates_smooth_means_stds_stes(loop,:,2);
%     diff_magdiffs_success_failure_mean_std_timenorm_rates(loop,:,3)=(success_timenorm_rates_smooth_means_stds_stes(loop,:,1)-failure_timenorm_rates_smooth_means_stds_stes(loop,:,1))./(success_timenorm_rates_smooth_means_stds_stes(loop,:,1)+failure_timenorm_rates_smooth_means_stds_stes(loop,:,1));
%     diff_magdiffs_success_failure_mean_std_timenorm_rates(loop,:,4)=(success_timenorm_rates_smooth_means_stds_stes(loop,:,2)-failure_timenorm_rates_smooth_means_stds_stes(loop,:,2))./(success_timenorm_rates_smooth_means_stds_stes(loop,:,2)+failure_timenorm_rates_smooth_means_stds_stes(loop,:,2));
% end    
% 
% absdiff_absmagdiffs_succ_fail_tnorm_crosscell_means_stds=zeros(54,4,2);
% 
% for i=1:54
%     for j=1:4
%         absdiff_absmagdiffs_succ_fail_tnorm_crosscell_means_stds(i,j,1)=mean(abs(diff_magdiffs_success_failure_mean_std_timenorm_rates(:,i,j)));
%         absdiff_absmagdiffs_succ_fail_tnorm_crosscell_means_stds(i,j,2)=std(abs(diff_magdiffs_success_failure_mean_std_timenorm_rates(:,i,j)));
%     end
% end
% 
% velocity correlations

% vel_by_rate=zeros(no_cells,10);
% for loop=1:no_cells
% tfile=load([datapath tfile_list_char(loop,:)]);
% cell_spikes_velcount=zeros(10,2);
% for j=1:10
%     x=find(vels_filt>((j-1)*4) & vels_filt<(((j-1)*4)+4));
%     for k=1:length(x)
%         xx=find(tfile>(pos(x(k),1)-.05) & tfile<(pos(x(k),1)+0.049));
%         cell_spikes_velcount(j,1)=length(xx);
%         cell_spikes_velcount(j,2)=cell_spikes_velcount(j,2)+1;
%     end
% end
% vel_by_rate(loop,:)=cell_spikes_velcount(:,1)./cell_spikes_velcount(:,2);
% end
% vel_by_rate(:,:)=vel_by_rate(:,:).*10;


%making a time-normalized 'velocity_map' for the time preceding starts by 2 sec
%to the time following ends by 2 sec

% success_timenorm_vels=zeros(length(success_markers_times),54);
% failure_timenorm_vels=zeros(length(failure_markers_times),54);
% 
% 
%     for i=1:length(success_markers_times)
%         for j=-10:-1
%             x=find(pos(:,1)>(success_markers_times(i,2)+(j*.2)) & pos(:,1)<(success_markers_times(i,2)+(j*.2)+.19999));
%             success_timenorm_vels(i,j+11)=mean(vels_filt(x));
%         end
%         for j=1:10
%             x=find(pos(:,1)>(success_markers_times(i,4)+((j-1)*.2)) & pos(:,1)<(success_markers_times(i,4)+(((j-1)*.2)+.19999)));
%             success_timenorm_vels(i,j+30)=mean(vels_filt(x));
%         end    
%         w=success_markers_times(:,3)-success_markers_times(:,2);
%         ww=success_markers_times(:,4)-success_markers_times(:,3);
%         trial_timebin=w(i)/10;
%         for j=1:10
%             x=find(pos(:,1)>success_markers_times(i,2)+((j-1)*trial_timebin) & pos(:,1)<(success_markers_times(i,2)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
%             success_timenorm_vels(i,j+10)=mean(vels_filt(x));
%         end
%         trial_timebin=ww(i)/10;
%         for j=1:10
%             x=find(pos(:,1)>success_markers_times(i,3)+((j-1)*trial_timebin) & pos(:,1)<(success_markers_times(i,3)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
%             success_timenorm_vels(i,j+20)=mean(vels_filt(x));
%         end
%     end
% 
%     for i=1:length(failure_markers_times)
%         for j=-10:-1
%             x=find(pos(:,1)>(failure_markers_times(i,2)+(j*.2)) & pos(:,1)<(failure_markers_times(i,2)+(j*.2)+.19999));
%             failure_timenorm_vels(i,j+11)=mean(vels_filt(x));
%         end
%         for j=1:10
%             x=find(pos(:,1)>(failure_markers_times(i,4)+((j-1)*.2)) & pos(:,1)<(failure_markers_times(i,4)+(((j-1)*.2)+.19999)));
%             failure_timenorm_vels(i,j+30)=mean(vels_filt(x));
%         end    
%         w=failure_markers_times(:,3)-failure_markers_times(:,2);
%         ww=failure_markers_times(:,4)-failure_markers_times(:,3);
%         trial_timebin=w(i)/10;
%         for j=1:10
%             x=find(pos(:,1)>failure_markers_times(i,2)+((j-1)*trial_timebin) & pos(:,1)<(failure_markers_times(i,2)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
%             failure_timenorm_vels(i,j+10)=mean(vels_filt(x));
%         end
%         trial_timebin=ww(i)/10;
%         for j=1:10
%             x=find(pos(:,1)>failure_markers_times(i,3)+((j-1)*trial_timebin) & pos(:,1)<(failure_markers_times(i,3)+((j-1)*trial_timebin)+(trial_timebin-.00000001)));
%             failure_timenorm_vels(i,j+20)=mean(vels_filt(x));
%         end
%     end
% 
% success_timenorm_vels_smooth=zeros(length(success_markers_times),54);
% failure_timenorm_vels_smooth=zeros(length(failure_markers_times),54);
% 
%     for i=1:no_successes
%         success_timenorm_vels_smooth(i,:)=filtfilt(y3_peth,1,success_timenorm_vels(i,:));
%     end
%     for i=1:no_failures
%         failure_timenorm_vels_smooth(i,:)=filtfilt(y3_peth,1,failure_timenorm_vels(i,:));
%     end

%take means stds stes of smoothed timenorm rates

% success_timenorm_vels_smooth_means_stds_stes=zeros(54,3);
% failure_timenorm_vels_smooth_means_stds_stes=zeros(54,3);
% 
%     for i=1:54
%     success_timenorm_vels_smooth_means_stds_stes(i,1)=mean(success_timenorm_vels_smooth(:,i));
%     success_timenorm_vels_smooth_means_stds_stes(i,2)=std(success_timenorm_vels_smooth(:,i));
%     success_timenorm_vels_smooth_means_stds_stes(i,3)=std(success_timenorm_vels_smooth(:,i))/(sqrt(no_successes-1));
%     failure_timenorm_vels_smooth_means_stds_stes(i,1)=mean(failure_timenorm_vels_smooth(:,i));
%     failure_timenorm_vels_smooth_means_stds_stes(i,2)=std(failure_timenorm_vels_smooth(:,i));
%     failure_timenorm_vels_smooth_means_stds_stes(i,3)=std(failure_timenorm_rates_smooth(:,i))/(sqrt(no_failures-1));
%     end

%take diffs and mag of diffs between timenorm success and failure patterns

% diff_magdiffs_success_failure_mean_std_timenorm_vels=zeros(54,4);
% 
%     diff_magdiffs_success_failure_mean_std_timenorm_vels(:,1)=success_timenorm_vels_smooth_means_stds_stes(:,1)-failure_timenorm_vels_smooth_means_stds_stes(:,1);
%     diff_magdiffs_success_failure_mean_std_timenorm_vels(:,2)=success_timenorm_vels_smooth_means_stds_stes(:,2)-failure_timenorm_vels_smooth_means_stds_stes(:,2);
%     diff_magdiffs_success_failure_mean_std_timenorm_vels(:,3)=(success_timenorm_vels_smooth_means_stds_stes(:,1)-failure_timenorm_vels_smooth_means_stds_stes(:,1))./(success_timenorm_vels_smooth_means_stds_stes(:,1)+failure_timenorm_vels_smooth_means_stds_stes(:,1));
%     diff_magdiffs_success_failure_mean_std_timenorm_vels(:,4)=(success_timenorm_vels_smooth_means_stds_stes(:,2)-failure_timenorm_vels_smooth_means_stds_stes(:,2))./(success_timenorm_vels_smooth_means_stds_stes(:,2)+failure_timenorm_vels_smooth_means_stds_stes(:,2));
%     
% 



clear  w* x* y* z* tfile trial_timebin tmp* pos_* rat_* m mm i j loop h hh hhh g gg ggg b* done ans 


        









