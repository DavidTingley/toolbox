function [clusters,order] = MV_spec(views,numclusts,graph_type,comb_type,clust_algo,plot_results,crossval)
% INPUTS
% views - Cell array (1xN) where N is the number of views to be used.  The
%       contents of each element in 'views' must be in the matrix form (M x D) 
%       where M equals the number of observations.  M must be the same values
%       across all elements in 'views'
%
% numclusts - Number of clusters (k-value) used when clustering
%
% graph_type - ('knn','gaussian','epsilon','cosine','corr') Variable determines how 
%              similarity graphs of each view are constructed.
%
% comb_type - ('sgt','ctn_comb','ctn_late') Variable determines how similarity graphs of each view are
%       combined. 
%
%       SGT - Joachims, T. (2003). Transductive learning via spectral graph partitioning. In Proceedings of the 20th inter-
%             national conference on machine learning (ICML 2003) (pp. 290–297).
%       sMD(ctn) - de Sa 2009
%       
%
% clust_algo- ('kmeans','hl','gmm') Variable determines which clustering algorithm will be used.
%       Options are K-means, Heirarchical Linkage, or a Gaussian Mixture Model
%
% plot_results - 0 to skip plotting results, 1 to view plotted results.
%
% crossval - 1 to skip cross validation, or any larger integer to choose the best
%       result out of that many iterations
%
% OUTPUTS 
% clusters - a vector of integer values that signifiy which category a
%       given observation is in
%
% TODO:
% - add cross validation
% - add documentation for sgt, ctn, ctn_late 
% - add JOINT(element wise mult) and matrix wise mult to graph combination
% - correlation or mahal dist for similarity graph


if ~iscell(views)
    disp('Warning: Views not in cell format, assuming single view input...')
    if isnumeric(views)  
    v{1} = views; clear views
    views = v; clear v
    end
end

if nargin < 2
	error('Must have at least the first two arguments:  MV_spec(views,numclusts)');
elseif nargin == 2
            graph_type = 'gaussian';  
            comb_type = 'sgt'; 
            clust_algo = 'kmeans';
            plot_results = 1;
            crossval = 1;
elseif nargin == 3
            comb_type = 'sgt'; 
            clust_algo = 'kmeans';
            plot_results = 1;
            crossval = 1;
elseif nargin == 4
            clust_algo = 'kmeans';
            plot_results = 1;
            crossval = 1;
elseif nargin == 5
            plot_results = 1;
            crossval = 1;
elseif nargin == 6
            crossval = 1;
end

%% check views
for t = 1:length(views)
    temp(t) = size(views{t},1);
end
if length(unique(temp)) > 1
    error('number of observations(rows) not consistant across views')
end

no_views=length(views);
randsweeps = 1000;
goto=size(views{1},1);

%% Graph Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch graph_type
    case 'knn'
        % KNN
        disp('Using KNN to make similarity graphs...')
        for j = 1:no_views
            A{j}=make_graph(views{j},'knn',round(size(views{1},1)*.1));
        end
        
    case 'gaussian'
        % Gaussian
        disp('Using a guassian kernel to make similarity graphs...')
        for j = 1:no_views
            A{j}=make_graph(views{j},'gaussianKernel');
        end
        
    case 'cosine'
        % cosine
        disp('Using a cosine similarity to make similarity graphs...')
        for j = 1:no_views
            A{j}=make_graph(views{j},'cosine');
        end
        
    case 'epsilon'
        % epsilon needs to be developed
        disp('Using an epsilon ball to make similarity graphs...')
        for j = 1:no_views
            A{j}=make_graph(views{j},'epsilonBall',round(size(views{1},1)*.1));
        end
        
    case 'corr'
        disp('Using correlation to make similarity graphs...')
        for j=1:no_views
           A{j} = corr(views{j},'rows','complete'); 
        end
end

for j = 1:no_views
    Azd{j} = A{j}-diag(diag(A{j}));
end


%% Graph Combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can this deal with missing data?-Groppe
switch comb_type
    case 'sgt'
        % sgt - Sums the similarity matrices produced from different views
        disp('Combining graphs using SGT...')
        SGT = zeros(size(views{1},1)); 
        
        for i = 1:no_views
            SGT = SGT + Azd{i};
        end
        
        DSGT=diag(sum(SGT'));
        LSGT=DSGT^(-.5)*SGT*DSGT^(-.5);
        [XSGTsort,DSGTsort]=eigs(LSGT,numclusts);
        XSGTuse=(XSGTsort(:,1:numclusts));
        XSGTsq=XSGTuse.*XSGTuse;
        divmatSGT=repmat(sqrt(sum(XSGTsq')'),1,numclusts);
        YSGT=XSGTuse./divmatSGT;
        
    case 'ctn_comb'
        % ctn normal
        disp('Combining graphs using CTN...')
        CTN_pairwise = zeros(size(views{1},1));
        
        for i = 1:no_views
            for j = i+1:no_views
                CTN_pairwise = A{i}*A{j} + CTN_pairwise;
            end
        end
        
        DCTN_pairbot = diag(sum(CTN_pairwise));
        DCTN_pairtop = diag(sum(CTN_pairwise'));
        LCTN_pair = DCTN_pairtop^(-.5)*CTN_pairwise*DCTN_pairbot^(-.5);
        [UCTN,SCTN,VCTN]=svd(LCTN_pair);
        XCTNcombuse=.5*UCTN(:,1:numclusts)+.5*VCTN(:,1:numclusts);
        XCTNcombsq=XCTNcombuse.*XCTNcombuse;
        divmatCTNcomb=repmat(sqrt(sum(XCTNcombsq')'),1,numclusts);
        YCTNcomb=XCTNcombuse./divmatCTNcomb;
        
    case 'ctn_late'
        % ctn late
        disp('Combining graphs with using CTN late...')
        CTN_pairwise = zeros(size(views{1},1));
        
        for i = 1:no_views
            for j = i+1:no_views
                CTN_pairwise = A{i}*A{j} + CTN_pairwise;
            end
        end
        
        DCTN_pairbot = diag(sum(CTN_pairwise));
        DCTN_pairtop = diag(sum(CTN_pairwise'));
        LCTN_pair = DCTN_pairtop^(-.5)*CTN_pairwise*DCTN_pairbot^(-.5);
        [UCTN,SCTN,VCTN]=svd(LCTN_pair);
        XCTNuse=[UCTN(:,1:numclusts);VCTN(:,1:numclusts)];
        XCTNsq=XCTNuse.*XCTNuse;
        divmatCTN=repmat(sqrt(sum(XCTNsq')'),1,numclusts);
        YCTN=XCTNuse./divmatCTN;
        YCTNlate=.5*YCTN(1:goto,:)+.5*YCTN(goto+1:goto*2,:);
end

%% Unsupervised Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Running Clustering Algorithm: ' clust_algo])

switch comb_type
    case 'ctn_comb'
        Y = YCTNcomb;
    case 'ctn_late'
        Y = YCTNlate;
    case 'sgt'
        Y = YSGT;
end

switch clust_algo
    case 'kmeans'
        [clusters,cls,distortion]=kmeans(Y,numclusts,'replicates',randsweeps);
    case 'hl'
        Z = linkage(Y,'single','mahalanobis');
        [clusters] = cluster(Z,'MaxClusts',numclusts);
    case 'gmm'
        [clusters P lo] = GMM(Y,numclusts);
        for i = 1:crossval-1
            [c,p,l] = GMM(Y,numclusts);
            if ~isnan(l)
                if l > lo | i == 1
                    clusters = c; P = p; L = l;
                    lo = l;
                end
            end
        end
end


%% Plotting
if plot_results == 1    
        clear temp
    for i = 1:numclusts
        temp(i,:) = mean(views{1}(clusters == i,:));   % I'm guessing views{3} should be views{1}-Groppe
        %ORIG LINE: temp(i,:) = mean(views{3}(clusters == i,:));
        [t t order] = sort_cells(temp,temp,1); %You need to add sort_cells.m to the toolbox--Groppe
    end
    

    
    p = factor(numclusts);
    m = p(end)+1;
    if length(p)>2
        n=0;
        for i = 1:length(p)-2
            n = n+times(p(i),p(i+1));
        end
    elseif length(p)==2
        n = p(end-1);
    elseif length(p)==1
        n = 1;
    end
    for f = 1:length(views)
        figure
        for i = 1:numclusts
            subplot(n,m,i)
            if size(views{f},2) > 1
                imagesc(views{f}(clusters == order(i),:))
            elseif size(views{f},2) == 1
                bar(views{f}(clusters == order(i),:))
            end
        end
    end
end

