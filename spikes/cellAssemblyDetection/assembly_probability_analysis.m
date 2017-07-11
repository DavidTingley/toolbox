


c = 1;
for i = 1:length(rec_place1)
% find B cell first
B = rec_place1(i);

AC = find(rec_place1==B);

for j = 1:length(AC)
for jj = 1:length(AC)
    if j ~= jj
    A = find(rec_place1==rec_place2(AC(j)));
    C = find(rec_place1==rec_place2(AC(jj)));
    allcombs = length(A)+length(C);
    uniquecombs = length(unique([A;C]));
    percentUnique(c) = uniquecombs/allcombs;
    
    ACControl = AC(randperm(length(AC)));
    CControl = find(rec_place2==rec_place1(ACControl(jj)));
    allcombsControl = length(A)+length(CControl);
    uniquecombsControl = length(unique([A;CControl]));
    percentUniqueControl(c) = uniquecombsControl/allcombsControl;
    
    c = c+1;
    end
end
end
end


