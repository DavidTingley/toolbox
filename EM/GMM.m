function [clusters,P,L,tries] = GMM(X,K,prin_comps)

%       This code executes an expectation maximization algorithm for 
%   multidimensional gaussian mixture model.  
%
%   INPUTS:
% X - (n x d) data matrix where rows(n) are observations
%
% K - integer that is the desired number of clusters
%
% prin_comps - integer less than d that determines the number of principle
%           components use.  This algorithm uses PCA to compress the data
%           set prior to clustering.
%
%   OUTPUTS:
% clusters - vector of values within 1:K signifying which cluster a given
%           observation belongs to
% P - (n x K) matrix of probabilities that a given observation is in a given
%           category

if nargin < 2
	error('Must have at least the first two arguments:  GMM(data,K)');
elseif nargin == 2
    prin_comps = round(size(X,2)/4);
    if prin_comps < 5
       prin_comps = size(X,2); 
    end
elseif nargin == 3

elseif nargin == 4

elseif nargin == 5

end

max_iterations = 1000;       % Max number of iterations
max_redos = prin_comps;
tries = 0;
t = cputime;
THRESHOLD = .10;  % or  change in log likelihood
num_loops = 0;
exit = 0;


[init_clust, init_place, score, W, Mus, Sigma] = EM_init(prin_comps, K,X);  
                                    %initializes Gaussian mean/covariances

[n d] = size(score(:,1:prin_comps)); 
L = EM_likelihood(score,W,Mus,Sigma);
Lo = 2*L;
warning off 

while exit == 0 

    %% E step = Compute responsibilities for each Gaussian for each waveform
    
   P = EM_EStep(K,d,n,W, Mus,Sigma,score);

    %% M step = Modify parameters of Gaussians to maximize likelihood for
    %          the weighted Data
     
    [W,Mus,Sigma] = EM_MStep(n,d,K,P,score);

   %%  Find new likelihood and check for exit parameters
    Lo = L;
    L = EM_likelihood(score,W,Mus,Sigma);
    num_loops    
    if num_loops >= max_iterations
        exit = 1;
        disp('Max Iterations reached...')
    end 
    if abs(100*(L-Lo)/Lo)<THRESHOLD 
       exit = 2; 
       disp('Converged...')
    end
    if isnan(L)
        if tries < max_redos      
        disp('Unable to converge: L = NaN, retrying with one less principle components...')
        tries = tries+1;
        [init_clust, init_place, score, W, Mus, Sigma] = EM_init(prin_comps-tries, K,X); 
        [n d] = size(score(:,1:prin_comps-tries)); 
        L = EM_likelihood(score,W,Mus,Sigma);
        Lo = 2*L;
        end
        if tries == max_redos
        disp('Unable to converge: L = NaN (try changing K or # of principle components)')
        disp('Using K-means clustering')
        clusters = init_clust;
        return;
        end
    end
    num_loops = num_loops + 1;
    
end

clusters = zeros(size(X,1),1);
for i = 1:size(P,1) 
   [d ii] = max(P(i,:)); 
   clusters(i) = ii;
end


% %% Plotting
% if exist('sort_cells.m','file')
%         clear temp
%     for i = 1:K
%         temp(i,:) = mean(X(clusters == i,:));   % I'm guessing views{3} should be views{1}-Groppe
%         %ORIG LINE: temp(i,:) = mean(views{3}(clusters == i,:));
%         [t t order] = sort_cells(temp,temp,1); %You need to add sort_cells.m to the toolbox--Groppe
%     end
%     
% 
%     
%     p = factor(K);
%     m = p(end)+1;
%     if length(p)>2
%         n=0;
%         for i = 1:length(p)-2
%             n = n+times(p(i),p(i+1));
%         end
%     elseif length(p)==2
%         n = p(end-1);
%     elseif length(p)==1
%         n = 1;
%     end
%   
%         figure
%         for i = 1:K
%             subplot(n,m,i)
%             if size(X,2) > 1
%                 [s s o] = sort_cells(X(clusters == order(i),:),X(clusters == order(i),:),1);
%                 imagesc(s)
%    
% %             line([43 43],[0 1500],'linewidth',2,'color','k')
% %             line([13 13],[0 1500],'linewidth',2,'color','k')
% %             line([21 21],[0 1500],'linewidth',2,'color','k')
% %             line([37 37],[0 1500],'linewidth',2,'color','k')
% %             caxis([.2 1.05])
%             elseif size(views{f},2) == 1
%                 bar(X(clusters == order(i),:))
%             end
%         end
% end

% for i = 1:K
% clusters(P(:,i) > .95) = i;
% end
elapsed_time = sprintf('CPU time used for EM: %5.2fs',cputime-t)