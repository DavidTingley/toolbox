function [cellinfo] = bz_mergeCellInfo(cellinfo1,cellinfo2)

% first attempt at merging cellinfo structs



fnames1 = fieldnames(cellinfo1);
fnames2 = fieldnames(cellinfo2);

if length(fnames1) ~= length(fnames2)
    error('are you sure these cellinfo files are the same type?')
end

cellinfo = struct('region',getfield(cellinfo1,'region'),...
                  'sessionName',getfield(cellinfo1,'sessionName'));

for i=1:length(fnames1)
    if isstruct(getfield(cellinfo1,fnames{i})) | ischar(getfield(cellinfo1,fnames{i}))
        cellinfo = setfield(cellinfo,fnames{i},getfield(cellinfo1,fnames{i}));
    else
        if ~strcmp(fnames1{i},'region') &  ~strcmp(fnames1{i},'sessionName') 
            cellinfo = setfield(cellinfo,fnames1{i},[getfield(cellinfo1,fnames1{i}), getfield(cellinfo2,fnames1{i})]);
        end
    end
end