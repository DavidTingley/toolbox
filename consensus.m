function [c,N,d] = consensus(range,data,numIter,cluster_tech)
% this code implements a consensus clustering algorithm, where a set of
% data is repeatedly clustered with different numbers of clusters(K) and 
% the results are compared to deteremine an "optimal" number of clusters. 
% Here, optimal is defined as consistent clustering results across multiple
% iterations with a maximal number of clusters(K)

% INPUT
%   range - 1:X vector of values that represent the number of clusters to
%           be used (i.e. a vector [2,3,4] will cluster the data with a K-value of 
%           2, 3, or 4 and then compare the results
%   data - M x D matrix of data to be clustered.  rows are observations
%   numIter - integer value that is the number of iterations of clustering
%           to be performed on each range element. Default = 50
%   cluster_tech - string that determines which clustering technique to be
%           used. Can be 'kmeans', 'gmm'
%
% OUTPUT
%   c - 
%   N - integer value selected from range that is the optimal number of
%       clusters for the given data.
%   d - vectory of optimization values across all of range. (higher is
%       better)


% TODO
% -add MVspec and others to available clustering techniques
% -rename variables
% -refine plotting functionality


%% check for inputs
    if nargin < 2
        z = data_generator(10);
        numIter = 50;
    end
    if nargin == 2
        z = data;
        numIter=10;
    end
    if nargin >= 3
       z = data; 
    end
    if nargin ~= 4
       cluster_tech = 'kmeans'; 
    end
    
%% Begin clustering
    c = zeros(size(z,1),size(z,1),length(range));
    switch cluster_tech
        case 'gmm'
            for i = 1:numIter
                for j = 1:length(range)
                    [clusters{i,j},P{i,j},L(i,j)] = GMM(z,range(j));
                    for ii = 1:size(z,1)
                    for jj = 1:size(z,1)
                    if clusters{i,j}(ii) == clusters{i,j}(jj) 
                    c(ii,jj,j) = 1 + c(ii,jj,j);
                    end
                    end  
                    end
                end
            end
        case 'kmeans'
            for i = 1:numIter
                for j = 1:length(range)
                    [clusters{i,j},idx{i,j},sumd{i,j}] = kmeans(z,range(j),'emptyAction','singleton');
                    for ii = 1:size(z,1)
                    for jj = 1:size(z,1)
                    if clusters{i,j}(ii) == clusters{i,j}(jj) 
                    c(ii,jj,j) = 1 + c(ii,jj,j);
                    end
                    end 
                    end
                end
            end    
    end

%% Plotting to determine K
for j = 1:length(range)
    figure(1)
        subplot(4,3,j)
        [temp ii] = sort(c(:,:,j));
        imagesc(c(:,:,j)); title(range(j))
    figure(2)
        subplot(4,3,j)
        hist(reshape(c(:,:,j),size(z,1)*size(z,1),1),20);title(range(j))
    d(j) = length(find(c(:,:,j)>.9*numIter | c(:,:,j) ==0));
    figure(3); subplot(4,3,j)
    imagesc(sortrows(sortrows(sortrows(c(:,:,j))')))
     title(range(j))
%     clustergram(squeeze(c(:,:,j)),'Colormap','redbluecmap')
end

figure(4)
plot(range,d)
[m i] = max(d);

% num = input(['Please pick a number of clusters(' num2str(range(i)) ' are recommended): ']);
% if ~isempty(num) 
%     N =  num;
% elseif isempty(num)
    N = range(i);
% end


end