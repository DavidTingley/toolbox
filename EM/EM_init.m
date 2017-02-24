function [init_clust, init_place, score, W, Mus, Sigma] = EM_init(prin_comps, K, X)
% This function initializes values for the EM algorithm


[n,d] = size(X(:,1:prin_comps));

% find the principal components

% [rlvm, frvals, frvecs, score, mn, dv] = pca(X);  % fucking matlab changed
% the pca outputs....
[COEFF, score, LATENT, TSQUARED, EXPLAINED] = pca(X);

%initialize with Kmeans algorithm 

[init_clust init_place init_sumd] = kmeans(score(:,1:d),K,'EmptyAction','singleton','replicates',500);
while ~isempty(find(init_sumd == 0))
[init_clust init_place init_sumd] = kmeans(score(:,1:d),K,'EmptyAction','singleton','replicates',500);
disp('Kmeans cluster with one observation(press enter to repeat Kmeans)')
% k = waitforbuttonpress

end
% Initialize the Means
Mus = init_place(:,1:prin_comps);

% Initialize the Covariance Matrices
% NOTE: init_sumd is a measure of distance between centroid and datapoints,
% it can be used as an initializer for the cov matrices of MOG

Sigma = zeros(prin_comps,prin_comps,K);  % first element = 1 to restrict cov
                                         % matrices to the diagonal values

 [D,I] = min(EM_sqdist(Mus', score(:,1:prin_comps)'));  % the variable "I"
                                                        % will be the Mean 
                                                        % that is closest to
                                                        % each individual point
for i = 1:K
   b = find(I == i);
   if length(b) > d*2;
       S = cov(score(b,1:d));
   else 
       S = cov(score(b,1:d));
   end
   Sigma(:,:,i) = S;
end
                                                        
                                                        
% Initialize the weights 

W = zeros(K,1);
for j = 1:K
   W(j) = length(find(I==j))/n; 
end


end

